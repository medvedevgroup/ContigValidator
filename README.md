This is a validation pipeline for checking the accuracy of an assembler. 

It uses simple exact matching using a suffix tree and bwa mem and outputs the alignment percentage with a reference genome. 

It also checks how many kmers were found in the assembly as compared to the reference genome. 

To get the kmer counts, we use dsk. 

Suffix tree is from - https://github.com/kvh/Python-Suffix-Tree

We also use bwa mem to align and get the percentage. 

INSTALLATION

To use the program

SSH - 
git clone git@github.com:mayankpahadia1993/validationpipeline.git


HTML - 
git clone https://github.com/mayankpahadia1993/validationpipeline.git

USAGE 

cd validationpipeline
bash run.sh -r referenceFile.fa -s suffixTreeOutput.p -i InputFile.fa -a AlignmentResults

OPTIONS

echo "-r -s -i|-f -a are compulsory options. -i can be multiple"
echo "Use -r | --reference filename -> for passing the reference file"
echo "Use -s | --suffixtree filename  -> filename will store the suffixtree"
echo "Use -i | --input filename -> filename contains the input files to be checked"
echo "Use -a | --alignment filename -> filename is the file in which alignment results will be shown"
echo "Use -f | --file filename -> filename is a file which contains the input filenames"


