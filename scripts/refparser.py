#!/usr/bin/env python


class BedParse(object):
    """bed file parser"""
    def __init__(self, bed):
        

        self.bed    = bed 
        self.parsed = self.parse()
        self.names  = [name for name in self.parsed]
        
    def parse(self):
        """Parse information from bed file to self.parsed[name] = [chr, lowest coordinate, highest coordinate, strand]"""
        result = {}
        for l in open(self.bed):
            chromosome, lowest_cor, hiest_cor, name, blank, entry_or = l.split("\t")
            result[name] = [chromosome, int(lowest_cor), int(hiest_cor), entry_or.strip("\n")]
        return result

    def get_chr(self, name):
        """Takes entry name and return chromosome"""
        if self.parsed.get(name):
            return self.parsed.get(name)[0]
        return "{} not found..".format(name)
        
    def get_start(self, name):
        """Takes entry name and return lowest coordinate"""
        if self.parsed.get(name):
            return self.parsed.get(name)[1]
        return "{} not found..".format(name)
    
    def get_stop(self, name):
        """Takes entry name and return highest coordinate"""
        if self.parsed.get(name):
            return self.parsed.get(name)[2]
        return "{} not found..".format(name)

    def get_strand(self, name):
        """Takes entry name and return if entry is located on +/- strand"""
        if self.parsed.get(name):
            return self.parsed.get(name)[3]
        return "{} not found..".format(name)
            
    def get_length(self, name):
        """Takes entry name and return length of transcript"""
        if self.parsed.get(name):
            return self.parsed.get(name)[2]-self.parsed.get(name)[1]
        return "{} not found..".format(name)

    def get_region(self, chromosome, start, stop):
        """return the names of entries within a give region """
        in_region = []
        for entry in self.parsed:
            entry_chrom = self.parsed[entry][0] == chromosome 
            entry_start = self.parsed[entry][1] >= start and self.parsed[entry][1] <= stop # True if entry start within given region
            entry_stop  = self.parsed[entry][2] >= start and self.parsed[entry][2] <= stop # True if entry start within given region
            entry_loc   = entry_start or entry_stop # True if entry within given region
            if entry_chrom and entry_loc:
                in_region.append(entry)
        return in_region
        

class GFFParse(object):
    """ gff file parser. Parse information based on the IDs in fasta_id """

    def __init__(self, gff, fasta_id):
        self.gff      = gff
        self.features = []
        self.name_id  = {}
        self.parents  = {}
        self.parsed   = self.parse(fasta_id)
        self.names    = [name for name in self.parsed]

        
        
    def parse(self, fasta_id):
        """Parse information from gff file to self.parsed[ddbID] = 
        [chromosome, lowest coordinate, highest coordinate, name, description]"""
        parsed = {}
        with open(self.gff) as fin:
            identified = 0
            for line in fin:
                line = line.split("\t")
                if len(line) != 9: continue

                if "JH723952.1_7" in line[8]: 1/0
                if "=" in line[8]:
                    info = [x[x.find("=")+1:].strip("\n") for x in line[8].split(";")]
                    for meta in info:
                        if meta in fasta_id:
                            if meta not in parsed:
                                identified += 1
                                parsed[meta] = [line[0], int(line[3]), int(line[4]), line[6], meta, "n/a"]
                else:
                    name = line[8].strip("\n")
                    if name in fasta_id:
                        if name not in parsed:
                            identified += 1
                            parsed[name] = [line[0], int(line[3]), int(line[4]), line[6], name, "n/a"]
        return parsed

    def parse_info(self, info_column):
        """ Retrieve ID, Name and description """
        ddbID = "n/a"
        name  = "n/a"
        desc  = "n/a"
        for entry in info_column.split(";"):
            if "ID=" in entry:
                ddbID = entry[3:].strip("\n")
            elif "Name=" in entry:
                name = entry[5:].strip("\n")
            elif "description=" in entry:
                desc = entry[12:].strip("\n")
        if name != "n/a":
            if name in self.name_id:
                print("Warning! Feature name: {} exist multiple times in gff.".format(name))
                self.name_id[name].append(ddbID)
            self.name_id[name] = [ddbID]
        return ddbID, name, desc

    def get_id(self, name):
        """ Return corresponding ddbID to name. Returns input if no ddbID is found """
        if self.name_id.get(name):
            return self.name_id[name]
        return "{} not found..".format(name)

    def get_chr(self, ddbID):
        """Takes entry ddbID and return chromosome"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[0]
        return "{} not found..".format(ddbID)
        
    def get_start(self, ddbID):
        """Takes entry ddbID and return lowest coordinate"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[1]
        return "{} not found..".format(ddbID)
    
    def get_stop(self, ddbID):
        """Takes entry ddbID and return highest coordinate"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[2]
        return "{} not found..".format(ddbID)

    def get_strand(self, ddbID):
        """Takes entry ddbID and return if entry is located on +/- strand"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[3]
        return "{} not found..".format(ddbID)

    def get_name(self, ddbID):
        """ Returns name associated with gene if any, otherwise return input"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[4]
        return "{} not found..".format(ddbID)

    def get_desc(self, ddbID):
        """ Returns description for the give ddbID """
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[5]
        return "{} not found..".format(ddbID)
            
    def get_length(self, ddbID):
        """Takes entry ddbID and return length of transcript"""
        if self.parsed.get(ddbID):
            return self.parsed.get(ddbID)[2]-self.parsed.get(ddbID)[1]
        return "{} not found..".format(ddbID)

    def get_region(self, chromosome, start, stop):
        in_region = []
        for entry in self.parsed:
            entry_chrom = self.parsed[entry][0] == chromosome 
            entry_start = self.parsed[entry][1] >= start and self.parsed[entry][1] <= stop # True if entry start within given region
            entry_stop  = self.parsed[entry][2] >= start and self.parsed[entry][2] <= stop # True if entry start within given region
            entry_loc   = entry_start or entry_stop # True if entry within given region
            if entry_chrom and entry_loc:
                in_region.append(entry)
        return in_region

# class geneIDGFFParse(object):

#     def __init__(self, gff):

#         self.gff      = gff
#         self.name_id  = {}
#         self.parsed   = self.parse()
#         self.names    = [name for name in self.parsed]

        
        
#     def parse(self):
#         """Parse information from gff file to self.parsed[ddbID] = 
#         [chromosome, lowest coordinate, highest coordinate, name, description]"""
#         parsed = {}
#         with open(self.gff) as fin:
#             for l in fin:
#                 l = l.split("\t")
#                 if len(l) != 9: continue
#                 name = l[8]
#                 contig = l[0]
#                 start = int(l[3])
#                 stop = int(l[4])
#                 strand = l[6]
#                 if name not in parsed:
#                     parsed[name] = [contig, start, stop, strand, name, "n/a"]
#                 else:
#                     if stop > parsed[name][2]:
#                         parsed[name][2] = stop
#                     else:
#                         continue
                    
#         return parsed


#     def get_id(self, name):
#         """ Return corresponding ddbID to name. Returns input if no ddbID is found """
#         if self.name_id.get(name):
#             return self.name_id[name]
#         return "{} not found..".format(name)

#     def get_chr(self, ddbID):
#         """Takes entry ddbID and return chromosome"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[0]
#         return "{} not found..".format(ddbID)
        
#     def get_start(self, ddbID):
#         """Takes entry ddbID and return lowest coordinate"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[1]
#         return "{} not found..".format(ddbID)
    
#     def get_stop(self, ddbID):
#         """Takes entry ddbID and return highest coordinate"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[2]
#         return "{} not found..".format(ddbID)

#     def get_strand(self, ddbID):
#         """Takes entry ddbID and return if entry is located on +/- strand"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[3]
#         return "{} not found..".format(ddbID)

#     def get_name(self, ddbID):
#         """ Returns name associated with gene if any, otherwise return input"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[4]
#         return "{} not found..".format(ddbID)

#     def get_desc(self, ddbID):
#         """ Returns description for the give ddbID """
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[5]
#         return "{} not found..".format(ddbID)
            
#     def get_length(self, ddbID):
#         """Takes entry ddbID and return length of transcript"""
#         if self.parsed.get(ddbID):
#             return self.parsed.get(ddbID)[2]-self.parsed.get(ddbID)[1]
#         return "{} not found..".format(ddbID)

#     def get_region(self, chromosome, start, stop):
#         in_region = []
#         for entry in self.parsed:
#             entry_chrom = self.parsed[entry][0] == chromosome 
#             entry_start = self.parsed[entry][1] >= start and self.parsed[entry][1] <= stop # True if entry start within given region
#             entry_stop  = self.parsed[entry][2] >= start and self.parsed[entry][2] <= stop # True if entry start within given region
#             entry_loc   = entry_start or entry_stop # True if entry within given region
#             if entry_chrom and entry_loc:
#                 in_region.append(entry)
#         return in_region

class IDConvert(object):
    """ Specific id name converter based on input fasta files used in OrthoFinder search """

    def __init__(self, fasta):

        self.id_orthoID, self.orthoID_id = self.parse(fasta)

    def parse(self, fasta):
        id_orthoID = {}
        orthoID_id = {}
        with open(fasta) as fin:
            for line in fin:
                if line[0] == ">":
                    line = line.split()
                    orthoID = line[0][1:]
                    if "|" in orthoID:
                        ID = orthoID[:orthoID.find("|")]
                    else: ID = orthoID
                    if ID in id_orthoID: 1/0
                    id_orthoID[ID] = orthoID
                    orthoID_id[orthoID] = ID
        return id_orthoID, orthoID_id

    def get_orthoID(self, id):
        return self.id_orthoID.get(id, "not found..")

    def get_ID(self, id):
        return self.orthoID_id.get(id, "not found..")


class Ortho(object):

    def __init__(self, tsv):

        self.orthos = self.get_orthos(tsv)

    def get_orthos(self, tsv):
        #org_index = {"ddi": 2, "asu": 1, "dfa": 3, "ppa": 7, "dla": 5, "dfi": 4, "dpu": 6}
        index = {2: "ddi", 1: "asu", 3: "dfa", 7: "ppa", 5: "dla", 4: "dfi", 6: "dpu"}
        orthos = {x: {} for x in index.values()}
        with open(tsv) as fin:
            next(fin)
            for line in fin:
                line = line.split("\t")
                for i in range(1, 8):
                    org = index[i]
                    for gene in line[i].split(", "):
                        if gene not in orthos[org]:
                            gene = gene.strip("\n")
                            orthos[org][gene] = {x: [] for x in index.values()}
                        for g in range(1,8):
                            #if g ==i:continue
                            orthos[org][gene][index[g]] += [x.strip("\n") for x in line[g].split(", ")]
                            orthos[org][gene][index[g]] = list(set(orthos[org][gene][index[g]]))
                        
        return orthos




class FastaIndex(object):

    def __init__(self, fai):

        self.length = self.parse(fai)

    def parse(self, in_file):
        parsed = {}
        with open(in_file) as fin:
            for line in fin:
                line = line.split("\t")
                contig = line[0]
                length = int(line[1])
                parsed[contig] = length
        return parsed

    def get_length(self, contig):
        return self.length.get(contig, "{} not found..".format(contig))

# Parsing reference information for organisms. Additional organisms can be added

ddi_bed = BedParse("annotations/ddi_curated.bed")
dfi_bed = BedParse("annotations/dfi_curated.bed")

ddi_id = IDConvert("ortho/d_dis.fasta")
dfi_id = IDConvert("ortho/d_fir_trimmed2.fasta")
orthologues = Ortho("ortho/Orthogroups.tsv")

ddi_gff = GFFParse("annotations/D_discoideum_190503.gff", ddi_id.id_orthoID)
dfi_gff = GFFParse("annotations/dfir_ASM27748v1.gff", dfi_id.id_orthoID)

ddi_fai = FastaIndex("fasta/dicty_chromosomal.fa.fai")
dfi_fai = FastaIndex("fasta/GCA_000277485.1_ASM27748v1_genomic.fna.fai")



