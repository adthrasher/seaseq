#!/usr/bin/env cwl-runner
cwlVersion: v1.0
baseCommand: [ liquidator.pl ]
class: CommandLineTool
label: BAM liquidator v1 on bam file for all metagenes
doc: perl liquidator.pl -g <gff file> -b <bam file> [-o <outfile name>]

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var var_output_name = function() {
      if (inputs.outfile == ""){
        if (inputs.bamfile != null) { return inputs.bamfile.nameroot; }
      }
   };

inputs:
  bamfile:
    type: File
    secondaryFiles: $(self.basename+".bai")
    inputBinding:
      prefix: -b

  gfffile:
    type: File
    inputBinding:
      prefix: -g

  outfile:
    type: string?
    inputBinding:
      prefix: -o
      valueFrom: |
        ${
            if (self == ""){
              return var_output_name();
            } else {
              return self;
            }
        }
    default: ""

outputs:
  promoters:
    type: File
    outputBinding: 
      glob: '*-promoters.pdf'

  genebody:
    type: File
    outputBinding:
      glob: '*-entiregene.pdf'

  promotersheatmap:
    type: File
    outputBinding: 
      glob: '*heatmap.promoters.png'

  genebodyheatmap:
    type: File
    outputBinding:
      glob: '*heatmap.entiregene.png'