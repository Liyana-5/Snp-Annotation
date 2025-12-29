""" Script to annotate SNP variants from a VCF file"""
import argparse
from weakref import ref
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import logging
import os
import gffutils
import vcfpy 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Configure logging to log into a file
logger = logging.getLogger()
logger.setLevel(logging.INFO)
fh = logging.FileHandler('3043418.log')
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)

# Initiate A class for annotating SNP 
class SNPAnnotation:
    # Initiate pahs to the files used 
    def __init__(self, vcf_file, gff_file, fasta_file):
        self.vcf_file = vcf_file
        self.gff_file = gff_file
        self.fasta_file = fasta_file
        
    
    def load_fasta(self):
        self.genome_sequence= {}
        #Load the FASTA sequence into a dictionary called genome_sequence
        try:
            with open(self.fasta_file, "r") as fasta:
                for record in SeqIO.parse(fasta, "fasta"):
                    self.genome_sequence[record.id] = record.seq
        except Exception as e:
            logger.error(f"Error loading FASTA file: {e}")
            exit(1)
        return self.genome_sequence

    def create_gff_db(self):
        #Create a GFF database in memory using gffutils the logic is from gffutils API documentation
        # Created in memory since assesment documentation does not ask to outbut a .db file
        try:
            self.gff_db = gffutils.create_db(self.gff_file, dbfn=":memory:", force=True, keep_order=True)
            logger.info("GFF database created in memory.")
        except Exception as e:
            logger.error(f"Error creating GFF database: {e}")
            exit(1)
        return self.gff_db
        
    
    def create_bar_plot(self, output_file, plot_file="variant_types.png"):
        #Generate a bar plot for variant types (non-coding,non-synonymous and synonymous)
        try:
            # Load the output file containing variant data
            variant_data = pd.read_csv(output_file, sep="\t")
            
            # Create the bar plot
            sns.countplot(data=variant_data, x="Type", palette="pastel")
            plt.ylabel("Count")
            plt.title("Distribution of Variant Types")
            
            # Save the plot to a file
            plt.savefig(plot_file)
            plt.close()  
            logger.info(f"Bar plot saved to: {os.path.abspath(plot_file)}")
        except Exception as e:
            logger.error(f"Error generating bar plot: {e}")
        return plot_file
    
    # Process vcf data and save in an output file 
    def process_vcf(self, output_file,gff_db,genome_sequence):
        low_qual =0 
        results =  []
        try:
            try:
                vcf_reader = vcfpy.Reader.from_path(self.vcf_file)
            except Exception as e:
                logger.error(f"Error loading VCF file: {self.vcf_file}: {e}")
                exit(1)
            for record in vcf_reader: # Iterate through vcf file
                if record.QUAL is not None and record.QUAL <= 20:
                    low_qual += 1
                    continue # Skip low quality variants
                
                
                chrom, pos, ref = record.CHROM, record.POS, record.REF
                alts = [alt.value for alt in record.ALT]
                
                logger.debug(f"Processing variant at {chrom}:{pos} REF:{ref} ALT:{alts}")
                
                # Skip non-SNP variants (indels / multi-base alleles)
                if len(ref) != 1 or any(len(a) != 1 for a in alts):
                    logger.info(f"Skipping non-SNP variant at {chrom}:{pos} REF:{ref} ALT:{alts}")
                    
                    continue
                
                coding_region = False
                """query the coding regions which include position of variant (pos) within 
                their genomic coordinates since SNP has only one nucleotide start and end is pos"""
                
                for feature in self.gff_db.region(seqid=record.CHROM, start=record.POS, end=record.POS, featuretype="CDS"):
                    coding_region = True
                    parents = list(gff_db.parents(feature.id))
                    if not parents:
                        logger.warning(f"No parents found for CDS feature: {feature.id} (type: {feature.featuretype})")
                    for parent in parents:
                        transcript_id = parent.id #Get transcript ID of parents from gff file
                        if "." not in transcript_id: # Get only transcrips no genes 
                            continue
                        """ when i used featuretype as mRNA in my script i found around 138 variants are skipped from output
                        so i added a logger warning so i found some cds features did not have mrna as feature type i understood 
                        these are pesudogene transcripts so i asked chatgpt4.o how to include both pseudone transcripts
                        as well as mRNA transcripts in my output but only include transcripts in the output 
                        the logic to only query with feaure id and to skip gene id without "." is given by chatgpt """
                        
                        # Build complete coding sequence by adjusting the frame of cds and store start and end in a dict
                        coding_seq = ""
                        coding_regions = [] # Container to hold start and end positions of a coding segment
                        """" the frame determines the start of translation in a coding sequence the if frame = 0
                        the translation starts at first nuclotide of CDS if 1 at second nucleotid and if 2 at 3 rd
                        so start is adjusted according to frame information from gff file for both negative and positive
                        strand """
                        if parent.strand == "+":
                            for child in self.gff_db.children(parent.id, featuretype="CDS", order_by="start"):
                                phase = child.frame
                                start = child.start
                                end = child.end
                                if phase == 1:
                                    start +=1
                                elif phase == 2:
                                    start +=2
                                coding_seq += genome_sequence[child.seqid][start-1:end] # -1 used to account for 0 based indexing
                                coding_regions.append((start, end))
                        else:
                            for child in self.gff_db.children(parent.id, featuretype="CDS", order_by="start", reverse=True):
                                phase = child.frame
                                start = child.start
                                end = child.end
                                if phase == 1:
                                    start -=1
                                elif phase == 2:
                                    start -=2
                                coding_seq += genome_sequence[child.seqid][start-1:end]
                                coding_regions.append((start, end))
                            coding_seq = Seq(coding_seq).reverse_complement()
                        
                        # Find the position of the variant in the coding sequence
                        """ the variant position is how many bases into the sequence the pos lies"""
                        coding_seq_length = 0
                        variant_positions =[] # Dict to store variant positions for coding sequences 

                        if parent.strand == "+":
                            for start, end in coding_regions:
                                if start <= pos <= end:
                                    variant_position = coding_seq_length + (pos -start)
                                    variant_positions.append(variant_position)
                                coding_seq_length += (end -start+ 1) # Update the length
                        else:
                            for start, end in reversed(coding_regions): # Reversed since negative strand is read from end to start
                                if start <= pos <= end:
                                    variant_position = coding_seq_length + (end - pos)
                                    variant_positions.append(variant_position)
                                coding_seq_length += (end - start + 1)
                            
                        # Find the codon and mutate
                
                        """ A codon contains three nucleotides thus integer division  
                        of variant_position by 3 gives the index of the codon in sequence
                        and multiplying by 3 gives the codon start """
                        for variant_position in variant_positions:
                            codon_start = (variant_position // 3) * 3
                            ref_codon = coding_seq[codon_start:codon_start + 3] # Get effected codon
                            # To mutate the codon we convert ref codon to seq object and stores the alternate codon in mutable object
                            alt_codon_mutable = MutableSeq(Seq(ref_codon))
                            for alt in alts:
                                # Mutate the codon
                                """ Mutate at the index of variant position given by
                                modulus division of variant position bu 3 in codon since a codon contain 3 neucleotide"""
                                if parent.strand == "+":
                                    alt_codon_mutable[variant_position % 3] = alt
                                else:
                                    # Since coding seq is reverse complemented reverse complement alt as well to get correct amino acid 
                                    alt_codon_mutable[variant_position % 3] = str(Seq(alt).reverse_complement()) 
                                
                                # Change alternate codon back into seq objec to translate
                                alt_codon = Seq(alt_codon_mutable)

                                # Get reference and alternate amino acids
                                ref_aa = (ref_codon).translate()
                                alt_aa = (alt_codon).translate()

                                """AA position is codon index of the codon containing variant in the sequence 
                                since codon index is 0 based we add one""" 
                                aa_position = variant_position // 3 + 1  
                                # Append to results 
                                results.append([chrom, pos, ref, alt, 'Non-synonymous' if ref_aa != alt_aa else 'Synonymous',
                                        transcript_id, aa_position, ref_aa, alt_aa if ref_aa != alt_aa else 'NA'])
                                
                if not coding_region :
                    for alt in alts:
                        results.append([chrom, pos, ref, alt, 'Non-coding', 'NA', 'NA', 'NA', 'NA'])
            
            # Create data frame to store into out table 
            data_frame = pd.DataFrame(results, columns=["CHROM", "POS", "REF", "ALT", "Type", 
                                                        "Transcript", "Protein Location", "Ref AA", "Alt AA"])
            data_frame.to_csv(output_file, sep="\t", index=False)
            logger.info(f"Number of variants with QUAL <= 20: {low_qual}") # Log number of low quality variants
            logger.info(f"Output table saved to {os.path.abspath(output_file)}")# Log path to ouput file
        except Exception as e:
            logger.error(f"Error processing VCF file: {e}")
        return output_file

    # Implement the methods in order 
    def execute(self):
        try:
            genome_sequence = self.load_fasta()
            gff_db = self.create_gff_db()
            output_file = "3043418_assesment.tsv"
            self.process_vcf(output_file, gff_db, genome_sequence)
            self.create_bar_plot(output_file)
        except Exception as e:
            logger.error(f"Execution failed: {e}")
            exit(1)
            
def main():
    #Use argparse to obtain files from the command line
    parser = argparse.ArgumentParser(description="Process VCF, GFF, and FASTA files")
    parser.add_argument('vcf_file', help="Path to the VCF file")
    parser.add_argument('gff_file', help="Path to the GFF file")
    parser.add_argument('fasta_file', help="Path to the FASTA file")
    args = parser.parse_args()

    # Log input file names
    logger.info(f"VCF file: {args.vcf_file}")
    logger.info(f"GFF file: {args.gff_file}")
    logger.info(f"FASTA file: {args.fasta_file}")

    # Check if files exists
    for file_path in [args.vcf_file, args.gff_file, args.fasta_file]:
        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            exit(1)
    #call the class and execute the class 
    annotation = SNPAnnotation(args.vcf_file, args.gff_file, args.fasta_file)
    annotation.execute()
    
if __name__ == "__main__":
    main()
