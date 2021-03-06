#!/usr/bin/env python
# Eric Tran
# Motif Mark OOP Assignment
# Bi 625 Winter 2021

import re
import itertools
import cairo
import argparse
import copy

class Motif:
    def __init__(self, motif_name):
        self.motif_name = motif_name
        self.length = 0

class FastaHeader:
    def __init__(self, header):
        self.header = header
        self.sequence = ""

class Gene:
    def __init__(self):
        self.sequence = ""
        self.length = 0
        self.drawing_coordinates = {}

    def set_length(self):
        self.length = len(self.sequence)

class Exon:
    def __init__(self):
        self.sequence = ""
        self.length = 0
        self.coordinates = {}

    def set_length(self):
        self.length = len(self.sequence)

    def set_exon_coordinates(self, sequence):
        capital_char_indices = []
        exon_coordinates = {}

        for count, char in enumerate(sequence):
            if char.isupper():
                capital_char_indices.append(count)
        exon_coordinates["exon"] = (capital_char_indices[0], capital_char_indices[-1])

        self.coordinates = exon_coordinates

class GeneGroup:
    def __init__(self, geneObject, exonObject, motifObject, fastaHeaderObject):
        self.gene = geneObject
        self.exon = exonObject
        self.motif = motifObject
        self.fasta = fastaHeaderObject

class CairoContext:
    def __init__(self, geneGroupObjectList):
        self.gene_group = geneGroupObjectList
        self.motif_list = []
        self.longest_intron = 0
        self.coordinates_dictionary = {}
        self.sequence_dictionary = {}
    
    def set_attributes(self):
        for i in range(4):
            header = self.gene_group[i].fasta.header
            sequence = self.gene_group[i].gene.sequence
            coordinates = self.gene_group[i].gene.drawing_coordinates

            self.coordinates_dictionary[header] = coordinates
            self.sequence_dictionary[header] = sequence

            # Set motif list
            self.motif_list.append(self.gene_group[i].motif.motif_name)

        for group in self.gene_group:
            if group.gene.length > self.longest_intron:
                self.longest_intron = group.gene.length
    
    def draw(self):
        # Calculate the longest motif
        longest_motif = 0
        for motif in self.motif_list:
            length = len(motif)
            if length > longest_motif:
                longest_motif = length
        
        # Calculate Width and Height of the surface
        longest_motif_font = 12 * longest_motif
        left_border = longest_motif_font + 37
        gene_count = len(self.coordinates_dictionary)
        height = gene_count * 150
        width = left_border + 30 + self.longest_intron

        with cairo.SVGSurface("motif_visualization.svg", width, height) as surface:
            cxt = cairo.Context(surface)

            # Generate a motif color pallete.
            color_palette = {}
            color_palette['ygcy'] = [1, 0.2, 0.2, 0.6] # Peach
            color_palette['GCAUG'] = [4, 0, 4, 0.5] # Purple
            color_palette['catag'] = [0.0, 1, 0.0, 1] # Green
            color_palette['YYYYYYYYYY'] = [0.0, 0.0, 1, 1] # Blue

            motif_index = 0
            # Create a legend for the motifs and a box representing their colors.
            for name in self.motif_list:
                if name == 'GCATG':
                    continue
                if name == 'cauag':
                    continue

                motif_height = motif_index * 35
                cxt.set_source_rgba(color_palette[name][0], color_palette[name][1], color_palette[name][2], color_palette[name][3])
                cxt.set_line_width(3)
                cxt.rectangle(10, motif_height, 15, 15)
                cxt.fill()

                cxt.set_source_rgb(0.1, 0.1, 0.1)
                cxt.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                cxt.set_font_size(15)
                cxt.move_to(28, motif_height + 10)
                cxt.show_text(name)
                motif_index += 1

            gene_index = 0
            # Draw an image for each gene from the FASTA input file.
            for key in self.sequence_dictionary:

                # Draw intron line (black line)
                sequence_length = len(self.sequence_dictionary[key])
                gene_height = 100 + (gene_index * 150)
                cxt.set_source_rgb(0.1, 0.1, 0.1)
                cxt.set_line_width(3)
                cxt.move_to(left_border, gene_height)
                cxt.line_to(left_border + sequence_length, gene_height)
                cxt.stroke()

                # Mark the exon location
                exon_coordinates_tuple = self.coordinates_dictionary[key]["exon"]
                exon_start = exon_coordinates_tuple[0]
                exon_end = exon_coordinates_tuple[1]
                cxt.set_source_rgb(0.1, 0.1, 0.1)
                cxt.set_line_width(3)
                exon_width = exon_end - exon_start
                cxt.rectangle(left_border + exon_end, gene_height - 20, exon_width, 40)
                cxt.stroke()

                # Create a title for each gene using its header line.
                header_height = 45 + (gene_index * 150)
                cxt.set_source_rgb(0.1, 0.1, 0.1)
                cxt.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                cxt.set_font_size(25)
                cxt.move_to(left_border, header_height)
                cxt.show_text(key)

                # Mark the motifs.
                self.motif_coordinates = self.coordinates_dictionary[key]
                for motif in self.motif_coordinates:
                    if motif == "exon":
                        continue
                    if motif == 'GCATG':
                        motif = 'GCAUG'
                    if motif == 'cauag':
                        motif = 'catag'
                    for motif_tuple in self.motif_coordinates[motif]:
                        cxt.set_source_rgba(color_palette[motif][0], color_palette[motif][1], color_palette[motif][2], color_palette[motif][3])
                        cxt.set_line_width(3)
                        cxt.move_to(left_border + motif_tuple[0], gene_height - 20)
                        cxt.line_to(left_border + motif_tuple[0], gene_height  + 20)
                        cxt.stroke()
                gene_index += 1
    
    # Call the create_drawing function.
    # create_drawing(coordinates_dictionary, longest_intron, motif_list, sequence_dictionary)

def argument_parser():
    '''Returns parser object'''
    # Create parser object.
    parser = argparse.ArgumentParser()

    # Create flags for motif file and fasta file.
    parser.add_argument('-f', '--fasta', type=str, required=True, help='Absolute file path \
        for fasta file')
    parser.add_argument('-m', '--motif', type=str, required=True, help='Absolute file path \
        for motif file)')
    return parser.parse_args()

def get_motif_coordinates(motif, sequence):
    '''Returns a list of tuples specifying the motif coordinates'''
    coordinate_tuples = []

    for match in re.finditer(motif, sequence, re.IGNORECASE):
        current_motif_coordinate = (match.start(), match.end() - 1)
        coordinate_tuples.append(current_motif_coordinate)
    return coordinate_tuples

def get_exon_coordinates(sequence):
    '''Returns a dictionary of the exon coordinates'''
    capital_char_indices = []
    exon_coordinates = {}

    for count, char in enumerate(sequence):
        if char.isupper():
            capital_char_indices.append(count)
    exon_coordinates["exon"] = (capital_char_indices[0], capital_char_indices[-1])

    return exon_coordinates

def get_ambigious_coordinates(motif_list, sequence):
    '''Returns a dictionary of ambigious reads and their coordinates'''

    # A list of all ambigious motifs.
    ambiguous_motif_list = []
    pattern = re.compile('^[ATGCUatgcu]+$')

    for motif in motif_list:
        if re.search(pattern, motif):
            continue
        else:
            ambiguous_motif_list.append(motif)
    
    # Now generate a dictionary of ambigious motif combinations and their coordinates
    # Key: Ambigious Motif
    # Value: Dictionary (key: motif, value: list of coordinate tuples)
    vIUPAC = {
    "A":["A"            ],
    "C":[    "C"        ],
    "G":[        "G"    ],
    "T":[            "T"],
    "U":[            "U"],
    "W":["A",        "T"],
    "S":[    "C","G"    ],
    "M":["A","C"        ],
    "K":[        "G","T"],
    "R":["A",    "G",   ],
    "Y":[    "C",    "T"],
    "B":[    "C","G","T"],
    "D":["A",    "G","T"],
    "H":["A","C",    "T"],
    "V":["A","C","G",   ],
    "N":["A","C","G","T"],
    "Z":[               ],
    }

    ambiguous_motif_coordinates = {}

    # The following piece of code was referenced from stackoverflow and modified.
    for ambiguous_motif in ambiguous_motif_list:

        groups = itertools.groupby(ambiguous_motif.upper(), lambda char:char not in vIUPAC)
        splits = []
        for b,group in groups:
            if b:
                splits.extend([[g] for g in group])
            else:
                for nuc in group:
                    splits.append(vIUPAC[nuc])
        combinations = [''.join(p) for p in itertools.product(*splits)]

        temporary_motif_coordinates = {}

        for motif in combinations:
            temporary_motif_coordinates[motif] = get_motif_coordinates(motif, sequence)
        
        all_coordinates = []

        for coordinate_list in temporary_motif_coordinates.values():
            all_coordinates += coordinate_list
        
        ambiguous_motif_coordinates[ambiguous_motif] = all_coordinates

    return ambiguous_motif_coordinates

def get_all_coordinates(motif_list, sequence):
    '''Returns a dictionary of motifs as keys and a list of their corresponding 
    coordinates as tuples'''

    # Key: motif
    # Value: List of tuples. tuple(beginning index, end index)
    motif_coordinates = {}
    exon_coordinates = get_exon_coordinates(sequence)
    ambigious_coordinates = get_ambigious_coordinates(motif_list, sequence)

    for motif in motif_list:
        motif_coordinates[motif] = get_motif_coordinates(motif, sequence)
    
    # Merge motif and exon dictionaries into one
    all_coordinates = {**motif_coordinates, **exon_coordinates}

    # Add ambiguous coordinates.
    for ambig in ambigious_coordinates:
        all_coordinates[ambig] = ambigious_coordinates[ambig]

    return all_coordinates

def main():

    # Retrieve parser object for fasta and motif file paths.
    args = argument_parser()
    fasta_file_path = args.fasta
    motif_file_path = args.motif

    motif_object_list = []
    motif_file = open(motif_file_path, 'r')

    # Convert motifs into motif object and store into list.
    for motif in motif_file:
        motif_name = motif.strip()

        # Create motif object.
        motif_object = Motif(motif_name)

        # Store motif objects into list.
        motif_object_list.append(motif_object)

    # Run through the motif list again. If there is a U base, create a motif
    # of the DNA sequence.
    ut_motif_object_list = []

    for obj in motif_object_list:
        motif_name = copy.deepcopy(obj.motif_name)

        if 'u' in motif_name:
            new_motif_name = motif_name.replace('u', 't')
            motif_object = Motif(new_motif_name)
            ut_motif_object_list.append(motif_object)
        elif 't' in motif_name:
            new_motif_name = motif_name.replace('t', 'u')
            motif_object = Motif(new_motif_name)
            ut_motif_object_list.append(motif_object)
        elif 'U' in motif_name:
            new_motif_name = motif_name.replace('U', 'T')
            motif_object = Motif(new_motif_name)
            ut_motif_object_list.append(motif_object)
        elif 'T' in motif_name:
            new_motif_name = motif_name.replace('T', 'U')
            motif_object = Motif(new_motif_name)
            ut_motif_object_list.append(motif_object)

    # Update the original motif list
    motif_object_list += ut_motif_object_list

    ut_motif_object_list.clear()

    # Read in Fasta file information.
    fasta_file = open(fasta_file_path, 'r')

    # List to store fasta header objects
    fasta_object_list = []
    # List to store gene objects
    gene_object_list = []
    # List to store exon objects
    exon_object_list = []
    # List to store GeneGroup objects
    gene_group_list = []

    # Extract header lines and create fasta header and Gene objects.
    for line in fasta_file:
        line = line.strip()
        
        if line[0] == '>':
            fasta_object = FastaHeader(line)
            gene_object = Gene()

            fasta_object_list.append(fasta_object)
            gene_object_list.append(gene_object)
            continue
        gene_object_list[-1].sequence += line
    
    # Set gene object length.
    for gene in gene_object_list:
        gene.set_length()

    # Create a list of the motif names.
    motif_list = []
    for motif in motif_object_list:
        motif_list.append(motif.motif_name)

    # Set drawing coordinates for each gene.
    for gene in gene_object_list:
        gene.drawing_coordinates = get_all_coordinates(motif_list, gene.sequence)

    # Create the Exon for each gene (uppercase letters). Loop through the genes and set Exon coordinates.
    for gene in gene_object_list:
        gene_sequence = gene.sequence
        exon_object = Exon()
        exon_object.set_exon_coordinates(gene_sequence)
        exon_object_list.append(exon_object)

    # Package everything into GeneGroup Objects. Store into GeneGroup object list.
    for i in range(4):
        gene_object = gene_object_list[i]
        exon_object = exon_object_list[i]
        motif_object = motif_object_list[i]
        fasta_header_object = fasta_object_list[i]
        gene_group_object = GeneGroup(gene_object, exon_object, motif_object, fasta_header_object)
        gene_group_list.append(gene_group_object)

    # Draw image.
    drawingObject = CairoContext(gene_group_list)
    drawingObject.set_attributes()
    drawingObject.draw()
main()