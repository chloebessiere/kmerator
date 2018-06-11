#!/bin/julia

# AUTHOR : Sébastien RIQUIER, IRMB, Montpellier

@everywhere using ParallelDataTransfer
using ArgParse
@everywhere using FastaIO
@everywhere using ProgressMeter


#Parse argument
s = ArgParseSettings()
@add_arg_table s begin
    "--unannotated", "-u"
      action = :store_true
      help = "activated if the provided initial fasta file correspond to an annotation external from Ensembl. Otherwise, use ensembl fasta files !"

      "--stringent", "-s"
        action = :store_true
        help = "if you think a gene-specific tag must be unique but present in ALL KNOWN TRANSCRIPT, if false, the kmer is considered as gene-specific is present only in variant of the corresponding gene, independently of the number know variant containing the kmer"

      "--verbose", "-v"
        action = :store_true
        help = "If you want this script talk too much"

    "--genome", "-g"
        help = "the genome fasta (.fa) or index by jellyfish for kmer request"
        #type = String
        required = true
    "--transcriptome", "-t"
        help = "the transcriptome fasta (.fa) (FASTA ONLY) for kmer request and transcriptional variant informations"
        #type = String
        required = true
    "--level", "-l"
        help = "type 'gene', 'transcript' or 'chimera' to extract specific kmer at these different levels. 'chimera' option must be done with 'unannotated' option!."
        required = true
    "--output", "-o"
        help = "directory of output"
        default = "."
    "--length"
        help = "length required for the kmer generation"
        default = 31
        arg_type = Int
    "fasta_file" #argument, no option
        help = "the fasta input file"
        required = true

    "--threshold"
        help = "FOR GENE LEVEL ONLY : Minimum fraction of annotated transcripts containing this kmer to admit it (default 0.5)"
        default = 0.5
        arg_type = Float64
end
parsed_args = parse_args(ARGS, s)

# variable attribution

output = parsed_args["output"]
if ismatch(r"\/$", output)
  output = replace(output, r"/$", "")
end


println("output directory: $output")



verbose_option = parsed_args["verbose"]
unannotated_option = parsed_args["unannotated"]
stringent_option = parsed_args["stringent"]
if stringent_option == true
println("stringent: yes")
end
admission_threshold = parsed_args["threshold"]
genome = parsed_args["genome"]
transcriptome = parsed_args["transcriptome"]
println("transcriptome : $transcriptome")
println("genome : $genome")

kmer_length = parsed_args["length"]
println("length: $kmer_length")
nbthreads =  nworkers()
println("nb process = $nbthreads")
fastafile = parsed_args["fasta_file"]
println("input sequences files = $fastafile")

level = parsed_args["level"]
println("level of kmer specificity = $level")
if unannotated_option == true
println("unannotated: yes")
end

if ! isfile("$genome")
  error("genome not found")
end
if ! isfile("$transcriptome")
  error("transcriptome not found")
end
if ! isfile("$fastafile")
  error("input fasta file not found")
end

#verif of fasta ensembl annotation
if ismatch(r".*\.fa", transcriptome) == true
  println("input fasta (*.fa) file, continue")
else
  error("error provided transcriptome : not a fasta file") # necessary ? i don't think so!
end

# create dictionary of transcriptome fasta



dico_transcriptome = Dict()

FastaReader("$transcriptome") do fr
for (desc, seq) in fr
desc_array = split(desc)
      # println( desc_array[1])
gene_name = replace(desc_array[7], "gene_symbol:", "")
ensembl_transcript_name = desc_array[1]
ensembl_gene_name = desc_array[4]

dico_transcriptome["$gene_name-$ensembl_transcript_name"] = seq
end
end


#directories creation
if ! isdir("$output/sequences")
    mkdir("$output/sequences")
end
if ! isdir("$output/sequences/$kmer_length")
    mkdir("$output/sequences/$kmer_length")
end
if ! isdir("$output/tags")
    mkdir("$output/tags")
end
if ! isdir("$output/tags/$kmer_length")
    mkdir("$output/tags/$kmer_length")
end
if ! isdir("$output/jellyfish_indexes")
    mkdir("$output/jellyfish_indexes")
end
if ! isdir("$output/jellyfish_indexes/$kmer_length")
    mkdir("$output/jellyfish_indexes/$kmer_length")
end


#verif of fasta file
if ismatch(r".*\.fa", fastafile) == true
  println("input fasta (*.fa) file, continue")
else
  error("input : not a fasta file")
end

#verif of genome/transcriptome are fasta, do the jellyfish count in this case
if ismatch(r".*\.fa", genome) == true
  run(`jellyfish count -C -m $kmer_length -s 10000 -t $nbthreads -o $output/jellyfish_indexes/$kmer_length/$(replace(basename(genome), ".fa", ".jl")) $genome`)
  genome = "$output/jellyfish_indexes/$kmer_length/$(replace(basename(genome), ".fa", ".jl"))"
  println("genome kmer index output : $genome")
end
if ismatch(r".*\.fa", transcriptome) == true
  run(`jellyfish count -m $kmer_length -s 10000 -t $nbthreads -o $output/jellyfish_indexes/$kmer_length/$(replace(basename(transcriptome), ".fa", ".jl")) $transcriptome`)
  transcriptome = "$output/jellyfish_indexes/$kmer_length/$(replace(basename(transcriptome), ".fa", ".jl"))"
  println("transcriptome kmer index output : $transcriptome")
end


#split each sequence of fasta input file into individual fastas
# WARNING work only with ensembl fasta descriptions!
if unannotated_option == false
    FastaReader("$fastafile") do fr
    for (desc, seq) in fr
      desc_array = split(desc)
      #println( desc_array[1])
      gene_name = replace(desc_array[7], "gene_symbol:", "")
      ensembl_transcript_name = desc_array[1]
      ensembl_gene_name = desc_array[4]
      if length("$seq") >= kmer_length
        println("$gene_name-$ensembl_transcript_name : good length, continue")
        FastaWriter("$output/sequences/$kmer_length/$gene_name-$ensembl_transcript_name.fa") do fwsequence
        #for (desc2, seq2) in fr
        write(fwsequence, [">$gene_name-$ensembl_transcript_name", "$seq"])
        end

      else println("$gene_name-$ensembl_transcript_name : wrong length, ignored")
      end
    end
    end
  # if the input fasta is not a part of ensembl annotaion
  else
    FastaReader("$fastafile") do fr
    for (desc, seq) in fr
    gene_name = desc
    if length("$seq") >= kmer_length
      FastaWriter("$output/sequences/$kmer_length/$desc.fa") do fwsequence
      write(fwsequence, [">$desc", "$seq"])
    end
  else println("$desc ========> wrong length, ignored")
    end
    end
end
end
println(readdir("$output/sequences/$kmer_length/"))

#println(keys(dico_transcriptome))

@passobj 1 workers() dico_transcriptome
@passobj 1 workers() unannotated_option
@passobj 1 workers() genome
@passobj 1 workers() transcriptome
@passobj 1 workers() level
@passobj 1 workers() output
@passobj 1 workers() kmer_length
@passobj 1 workers() stringent_option
@passobj 1 workers() admission_threshold


#Create the function f for extraction of specifics kmers for one sequence fasta
#file.

@everywhere function f(splitted_fasta_files)
 
sleep(myid()*5)

println("$splitted_fasta_files")
#println(dico_transcriptome)

# manage the name of transcript/gene
  if unannotated_option == false # if annotated
    gene_name = split(splitted_fasta_files, "-")[1]
#    println("$gene_name")
    transcript_name = split(splitted_fasta_files, "-")[2]
#    test = filter((k,v) -> startswith(k, "$gene_name-"), dico_transcriptome)
#    println(test)
    nb_variants = length(filter((k,v) -> startswith(k, "$gene_name-"), dico_transcriptome))
  else # if unannotated
    gene_name = "$splitted_fasta_files"
    transcript_name = "$splitted_fasta_files"
  end
println("level : $level")
if level == "gene"
  tag_file = "$gene_name-$transcript_name-gene_specific.fa"
end
if level == "transcript"
  tag_file = "$gene_name-$transcript_name-transcript_specific.fa"
end
if level == "chimera"
  tag_file = "$gene_name-chimera_specific.fa"
end


fasta_array = Array([])



# take the transcript sequence for jellyfish query
sequence_fasta = "error if seen"
FastaReader("$output/sequences/$kmer_length/$splitted_fasta_files") do fr
for (desc, seq) in fr
    println("$desc")
  sequence_fasta = "$seq"
end
end
#println("$sequence_fasta")
#println("")

#global genome
  kmercounts_genome = readstring(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$genome"`)
  println("jellyfish query on genome.jf finished")
  kmercounts_transcriptome = readstring(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$transcriptome"`)

  println("jellyfish query on transcriptome.jf finished")



kmercounts_genome = split(kmercounts_genome, "\n")[1:end-1]
kmercounts_genome_dico = Dict()
for mer in kmercounts_genome
mer = split(mer)
counts = mer[2]
seq = mer[1]
kmercounts_genome_dico["$seq"] = counts
end

kmercounts_transcriptome_dico = Dict()
kmercounts_transcriptome = split(kmercounts_transcriptome, "\n")[1:end-1]
for mer in kmercounts_transcriptome
mer = split(mer)
counts = mer[2]
seq = mer[1]
kmercounts_transcriptome_dico["$seq"] = counts
end



i = 0

@showprogress 10 "k-mer processing... $gene_name" for mer in keys(kmercounts_transcriptome_dico)

#println("-------------------------------------------------------------------------------------------")

if !any(x->x == "$mer", keys(kmercounts_genome_dico))

    #println("$mer not defined : reversed ?")
  d = Dict("A"=>"T", "C"=> "G", "G"=>"C", "T"=>"A")
  new_mer = reverse(replace(mer, r"[ACGT]{1}", x->d[x]))  # "TGCA"
  #println("$(String(reverse_complement!(dna"$mer"))))")
  genome_count = kmercounts_genome_dico["$new_mer"]
  #println("genome count : $genome_count")
else

genome_count = kmercounts_genome_dico["$mer"]

end

transcriptome_count = kmercounts_transcriptome_dico["$mer"]
#println("transcriptome count : $transcriptome_count")
if level == "gene"
         #println("genome count : $genome_count")
          if genome_count == "1" || genome_count == "0"#if the kmer is present/unique or does not exist (splicing?) on the genome
            mer_regex = Regex(mer)
            transcripts_containing_this_kmer = keys(filter((k,v) -> ismatch(mer_regex, v), dico_transcriptome))
            #if nb_variants != length(transcripts_containing_this_kmer)
                #println("$mer : not equal : transcriptome count $transcriptome_count variant number $nb_variants & transcripts_with this kmer $transcripts_containing_this_kmer")
            #end
            #println(" number_variants = $nb_variants")
            #println(transcripts_containing_this_kmer)
            if stringent_option == true && float(transcriptome_count) == float(nb_variants) == float(length(transcripts_containing_this_kmer))
                #println("specific kmer found !")
              i = i+1
                tmp = length(transcripts_containing_this_kmer)
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)")
              push!(fasta_array,"$mer")
              #write(fw, [">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)", "$mer"])

          elseif stringent_option == false && all(x->startswith(x, "$gene_name-"), transcripts_containing_this_kmer) == true && float(transcriptome_count) > float(nb_variants*admission_threshold)

                #println("specific kmer found")
                i = i+1
                tmp = length(transcripts_containing_this_kmer)
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)" )
              push!(fasta_array, "$mer")
              #write(fw, [">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)", "$mer"])

            elseif unannotated_option == true && float(transcriptome_count) == float(0)

              i = i+1
              push!(fasta_array,">$gene_name-$transcript_name.kmer$i" )
              push!(fasta_array,"$mer")
              #write(fw, [">$gene_name-$transcript_name.kmer$i", "$mer"])

          end #end of fasta writing function
        
          #FastaWriter("$output/tags/$kmer_length/$tag_file") do fw
 #         writedlm("$output/tags/$kmer_length/$tag_file", transpose(fasta_array), "\n")
        end
  end



  #println("$genome_count")
if level == "transcript"
  #println("$genome_count")
  #println("$transcriptome_count")
  #println(parse(Int, genome_count) <= 1 )
  #println(float(transcriptome_count) == float(0))
  if unannotated_option == true && float(transcriptome_count) == float(0) && parse(Int, genome_count) <= 1
      i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")

  elseif float(transcriptome_count) == float(1) && parse(Int, genome_count) <= 1
    i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")
  end
 end

if level == "chimera" && unannotated_option == true
    #println("$genome_count")
    #println("$transcriptome_count")

    #println(parse(Int, genome_count) == 0 )
    #println(float(transcriptome_count) == float(0))
  if float(transcriptome_count) == float(0) && parse(Int, genome_count) == 0
    #println("YES")
    i = i+1
              push!(fasta_array,">$gene_name.kmer$i")
              push!(fasta_array,"$mer")
  end
  end


end
println("$level $gene_name : $i specific kmers found")

writedlm("$output/tags/$kmer_length/$tag_file", reshape(fasta_array,length(fasta_array)), "\n")
end #end of function
splitted_fasta_files = readdir("$output/sequences/$kmer_length/")

#pmap(f, splitted_fasta_files, unannotated_option, genome, transcriptome, level, output, kmer_length, stringent_option, admission_threshold)
pmap(f, splitted_fasta_files)

