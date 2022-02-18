import argparse
import sys
import glob
import pysam as pysamsies
import subprocess
from multiprocessing import Pool
from multiprocessing import cpu_count

PARSER = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

# vcf file and out-directory
PARSER.add_argument("-v", "--vcf", required=True,
                    help="VCF file with parental strains")
PARSER.add_argument("-o", "--outdir", required=True,
                    help="Directory where output files will be written")

PARSER.add_argument("-s1", "--strain1", required=True,
                    help="strain 1 in ASE")
PARSER.add_argument("-s2", "--strain2", required=True,
                    help="strain 2 in ASE")

PARSER.add_argument("-b", "--bams", required=True,
                    help="list bams")

PARSER.add_argument("-t", "--tosser", required=False,
                    help="create BAM with unmapped reads", action="store_true")

PARSER.add_argument("-over", "--coverage_over", required=False, default=1.50,
                    help="Maximum coverage depth as a multiple of genome-wide average")
PARSER.add_argument("-under", "--coverage_under", required=False, default=0.25,
                    help="Minimum coverage depth as a multiple of genome-wide average")

PARSER.add_argument("-qds", "--qds", required=False, default=2,
                    help="Mininum variant score")
PARSER.add_argument("-sor", "--sor", required=False, default=3,
                    help="Maximum strand bias score")
PARSER.add_argument("-mq", "--mps", required=False, default=40,
                    help="Minimum mean mapping quality score")

PARSER.add_argument("-minq", "--minq", required=False, default=40,
                    help="Minimum RNA read mapping quality")

PARSER.add_argument("-mask", "--masking_file", required=False, default=None,
                    help="File with genomic regions to mask"
                         "Chrom\tbeg\tend\n for masking")
PARSER.add_argument("-f", "--min_scaffold", required=False, default=1,
                    help="Minimum chromsome/scaffold length")

###############
N_CPU = cpu_count()

ARGIES = PARSER.parse_args()
ARGDICT = {}

ARGDICT["vcf"] = ARGIES.vcf
ARGDICT["outdir"] = ARGIES.outdir

ARGDICT["strain1"] = ARGIES.strain1
ARGDICT["strain2"] = ARGIES.strain2

ARGDICT["bams"] = ARGIES.bams.split(",")

THE_STUFF = set([ARGDICT["strain1"], ARGDICT["strain2"]])

ARGDICT["tosser"] = ARGIES.tosser

if ARGIES.masking_file:
    ARGDICT["masking_file"] = ARGIES.masking_file
ARGDICT["min_scaffold"] = int(ARGIES.min_scaffold)
ARGDICT["coverage_over"] = float(ARGIES.coverage_over)
ARGDICT["coverage_under"] = float(ARGIES.coverage_under)

ARGDICT["qds"] = float(ARGIES.qds)
ARGDICT["mps"] = float(ARGIES.mps)
ARGDICT["sor"] = float(ARGIES.sor)
ARGDICT["minq"] = float(ARGIES.minq)

##############

# COVERAGE
def coverage():
    """Calculates average coverage per strain/individual"""
    print("ITERATING OVER VCF FILE TO GET COVERAGE INFORMATION")
    cov = {}
    final_cov = {}
    the_right_stuff = []
    openvcf = open(ARGDICT["vcf"])
    for line in openvcf:
        if line[0] == "#":
            if line.split("=")[0] == "##contig":
                line = line.split("=")
                vcf_scaff = line[2].split(",")[0]
                vcf_scaff_lens = (line[3].split(">")[0])
                if int(vcf_scaff_lens) > ARGDICT["min_scaffold"]:
                    the_right_stuff.append(vcf_scaff)
            elif line[0:6] == "#CHROM":
                header_strains = line.strip().split("\t")[9:]
                for strain in header_strains:
                    if strain in THE_STUFF:
                        cov[strain] = {}
                        cov[strain]["cov"] = 0.0
                        cov[strain]["total"] = 0.0
        else:
            line = line.strip().split("\t")
            info = line[9:]
            if (line[0] in the_right_stuff
                    and len(line[3]) == 1
                    and max([len(i) for i in line[4].split(",")]) == 1):
                for ix_strain in enumerate(header_strains):
                    if "./" not in info[ix_strain[0]] and ix_strain[1] in THE_STUFF:
                        reads = info[ix_strain[0]].split(":")[1].split(",")
                        parent_cov = sum([int(i) for i in reads])
                        cov[ix_strain[1]]["cov"] = cov[ix_strain[1]]["cov"]  + parent_cov
                        cov[ix_strain[1]]["total"] = cov[ix_strain[1]]["total"] + 1
    openvcf.close()
    for strain in cov:
        avecov = cov[strain]["cov"]/cov[strain]["total"]
        print("SAMPLE: %s\tCOVERAGE: %s"%(strain, round(avecov, 2)))
        final_cov[strain] = avecov
    return(the_right_stuff, final_cov)


def error(message):
    """Prints error messages"""
    sys.exit(message)


def line_parser(line):
    """Parses VCF string to get mapping quality info"""
    linedict = {}
    line = line.split(";")
    for keyval in line:
        keyval = keyval.split("=")
        if len(keyval) > 1:
            values = keyval[1].split(",")
            if len(values) == 1:
                linedict[keyval[0]] = float(keyval[1])
            else:
                linedict[keyval[0]] = min([float(i) for i in values])
    return(linedict)


def masker(scaffs):
    """Specifies regions of genome to mask"""
    maskdict = {}
    if "masking_file" in ARGDICT:
        maskfile = open(ARGDICT["masking_file"])
        for line in maskfile:
            line = line.split("\t")
            try:
                int(line[1]) and int(line[2])
            except ValueError:
                error("CANNOT CONVERT MASKING COORDS TO INTEREGRS. "
                      "EXITING PROGRAM.")
            if line[0] not in maskdict:
                maskdict[line[0]] = set(range(int(line[1]), int(line[2])+1))
            maskdict[line[0]] = maskdict[line[0]] | set(range(int(line[1]), int(line[2])+1))
    else:
        for scaff in scaffs:
            maskdict[scaff] = set([-1])
    return(maskdict)


def get_vcftuple(chroms, final_cov):
    """Parses the VCF file, filters SNPs, and extracts relevant linermation"""
    print("PARSING VCF TO ANALYZE VARIANTS")
    the_allele_outdict = {}
    all_pos_dict = {}
    masking = masker(chroms)
    openvcf = open(ARGDICT["vcf"], "r") # feeds the file line by line (see below)
    # this is to create new files in the folder and to overwrite existing ones if already there.
    for line in openvcf: # let's go over each linne
        if line[0] == "#":
            if line[0:6] == "#CHROM":
                headersies = line.rstrip().split("\t")[9:]
                strain_1_pos = headersies.index(ARGDICT["strain1"])
                strain_2_pos = headersies.index(ARGDICT["strain2"])
        line = line.rstrip().split("\t")
        chrom = line[0]
        # need to import the line dict file here
        if (chrom in chroms
                and len(line[3]) == 1
                and max([len(i) for i in line[4].split(",")]) == 1):
            if  chrom not in the_allele_outdict:
                the_allele_outdict[chrom] = {}
                all_pos_dict[chrom] = set()
            pos = int(line[1])
            quals = line_parser(line[7])
            passy = None
            ref_alts = line[3].split(",") + line[4].split(",")
            #print(line)
            if len(ref_alts) == 2 and "*" not in ref_alts:
                #print(pos, quals)
                if ("QD" in quals and "MQ" in quals and "SOR" in quals):
                    if (pos not in masking[chrom]
                            and quals["QD"] >= ARGDICT["qds"]
                            and quals["MQ"] >= ARGDICT["mps"]
                            and quals["SOR"] < ARGDICT["sor"]):
                        passy = True
                if passy:
                    #print("pass")
                    check = []
                    strain_info = line[9:]
                    if "./" not in strain_info[strain_1_pos] + strain_info[strain_2_pos]:
                        strainy_1_info = strain_info[strain_1_pos].split(":")
                        strainy_2_info = strain_info[strain_2_pos].split(":")
                        strainy_1_info_genotypo = set(strainy_1_info[0].split("/"))
                        strainy_2_info_genotypo = set(strainy_2_info[0].split("/"))
                        #print(strainy_1_info_genotypo, strainy_2_info_genotypo)
                        if (len(strainy_1_info_genotypo) == 1 and len(strainy_2_info_genotypo) == 1
                               and strainy_1_info_genotypo != strainy_2_info_genotypo):
                            #print("PASS 2")
                            for strainy in [strain_1_pos, strain_2_pos]:
                                reads = strain_info[strainy].split(":")[1].split(",")
                                coveragie = sum([int(i) for i in reads])
                                current = headersies[strainy]
                                norm_cov = final_cov[current]
                                #print(current, coveragie, norm_cov)
                                if norm_cov*ARGDICT["coverage_under"] <= coveragie <= norm_cov*ARGDICT["coverage_over"]:
                                    check.append(current)
                        if len(check) == len(THE_STUFF):
                            strainy_1_allela = ref_alts[int(list(strainy_1_info_genotypo)[0])]
                            strainy_2_allela = ref_alts[int(list(strainy_2_info_genotypo)[0])]
                            all_pos_dict[chrom].add(pos-1)
                            the_allele_outdict[chrom][pos-1] = {}
                            the_allele_outdict[chrom][pos-1][strainy_1_allela] = ARGDICT["strain1"]
                            the_allele_outdict[chrom][pos-1][strainy_2_allela] = ARGDICT["strain2"]
    return(all_pos_dict, the_allele_outdict)


def origins(read):
    chrom = read.reference_name
    if read.is_proper_pair and read.mapping_quality >= ARGDICT["minq"]:
        orgs = set()
        sequence = read.query_sequence
        where_aligns = read.get_reference_positions(full_length=True)
        if chrom in VCF_DICTS[0]:
            matcha = set(where_aligns) & VCF_DICTS[0][chrom]
            if len(matcha) >= 1:
                for matchup in sorted(matcha):
                    the_allele = sequence[where_aligns.index(matchup)]
                    if the_allele in VCF_DICTS[1][chrom][matchup]:
                        origin = VCF_DICTS[1][chrom][matchup][the_allele]
                        orgs.add(origin)
                    else:
                        orgs.add("Unknown")
    else:
        orgs = set(["Fail"])
    return(orgs)

def names(read):
    return(read.query_name)

def process_reads(chrom_bam):
    chrom = chrom_bam["chrom"]
    bam = chrom_bam["bam"]
    bam_prep = bam.split("/")[-1].split(".")[0]
    bam_parse = pysamsies.AlignmentFile(bam, "rb")
    ###
    read_stuff = bam_parse.fetch(contig=chrom, start=0)
    originals = [origins(read) for read in read_stuff]
    read_stuff = bam_parse.fetch(contig=chrom, start=0)
    namies = [names(read) for read in read_stuff]
    ###
    reads_out = {}
    for read_info in zip(namies, originals):
        if read_info[0] not in reads_out:
            reads_out[read_info[0]] = read_info[1]
        else:
            reads_out[read_info[0]] = reads_out[read_info[0]] | read_info[1]
    ###
    if ARGDICT["tosser"]:
        processed_reads = set()
    for strain in THE_STUFF:
        outbam = pysamsies.AlignmentFile("{}/{}_{}_{}_temp.bam".format(ARGDICT["outdir"], bam_prep, strain, chrom), "wb", template=bam_parse)
        read_names = set([read_name for read_name in reads_out if reads_out[read_name] == set([strain])])
        read_stuff = bam_parse.fetch(contig=chrom, start=0)
        reads_to_output = set([read for read in read_stuff if read.query_name in read_names])
        for readsies in reads_to_output:
            outbam.write(readsies)
        if ARGDICT["tosser"]:
            processed_reads = processed_reads | read_names
    if ARGDICT["tosser"]:
        outbam = pysamsies.AlignmentFile("{}/{}_{}_{}_temp.bam".format(ARGDICT["outdir"], bam_prep, "Tosser", chrom), "wb", template=bam_parse)
        read_stuff = bam_parse.fetch(contig=chrom, start=0)
        unprocessed_reads = set([read for read in read_stuff if read.query_name not in processed_reads])
        for readsies in unprocessed_reads:
            outbam.write(readsies)

def process_bam(bam):
    print("processing {}".format(bam))
    bam_prep = bam.split("/")[-1].split(".")[0]
    chrom_bam = []
    for chrom in CHROMS:
        chrom_bam.append({"chrom": chrom, "bam": bam})
    # function keys here
    p = Pool(N_CPU)
    p.map(process_reads, chrom_bam)
    all_mergers = THE_STUFF | set (["Tosser"])
    for ds in all_mergers:
        temp_bams = glob.glob("{}/{}_{}*_temp.bam".format(ARGDICT["outdir"], bam_prep, ds))
        subprocess.run("samtools merge {}/{}_{}_merged.bam {}".format(ARGDICT["outdir"], bam_prep, ds, " ".join(temp_bams)), shell=True)
        subprocess.run("samtools sort -@ {} -T /tmp/{}.sorted -o {}/{}_{}.bam {}/{}_{}_merged.bam && samtools index {}/{}_{}.bam".format(N_CPU, ds, ARGDICT["outdir"], bam_prep, ds, ARGDICT["outdir"], 
                                                                                                                             bam_prep, ds, ARGDICT["outdir"], bam_prep, ds), shell=True)
    subprocess.run("rm {}/*_temp.bam && rm {}/*_merged.bam".format(ARGDICT["outdir"], ARGDICT["outdir"]), shell=True)

######################

STUFF = coverage()
CHROMS = STUFF[0]
FINAL_COV = STUFF[1]

VCF_DICTS = get_vcftuple(CHROMS, FINAL_COV)

for bammy in ARGDICT["bams"]:
    process_bam(bammy)

##################
