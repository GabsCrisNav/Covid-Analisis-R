#Librerias
library(Biostrings)
library(seqinr)
library(ade4)
library(ape)
library(ggtree)
library(bioseq)
library(DECIPHER)
library(ggplot2)
library(gridExtra)

#Funciones
DNA2b <- function(dna){
  
  dnal = strsplit(dna,"")
  dnab = as.DNAbin(dnal) 
  return(dnab)
  
}
#Convierte las secuencias de una archivo fasta en formato DNABin
nuc <- function(dna){
  
  l <- vector();
  
  for (i in 1:length(dna)) {
    seq <- dna[i]
    l<- c(l,count(seq[[1]], 1))
  }
  
  return(l)
}
#Obtiene las concentraciones de nucleotidos de una secuencia
Nams <- function(names){
  
  l <- vector()
  
  for (i in 1:length(names)) {
    gn <- names[i]
    l<- c(l, rep((c(gn)),c(4)))
  }
  
  return(l)
  
}
#Genera un vector con los nombre cuadriplicados de las secuencias

#Obtención de las secuencias

#Vector con indeficadores de los documentos Fasta
virus <- c(  "JX869059", "AY508724", "MN908947", "AY390556", "AY278489", "MN985325","AY485277","MT292571")

#Lee los atributos de las cadenas asignadas a los identificadores
virus_sequences <- read.GenBank(virus)

#Genera un Fasta que combina todas las secuencias
write.dna(virus_sequences,  file ="virus_seqs.fasta", format = "fasta", append =
            FALSE, nbcol = 6, colsep = " ", colw = 10)

#Visualiza el arreglo de las secuencias sin alinear con biostrings
virus_seq_not_align <- readDNAStringSet("virus_seqs.fasta", format = "fasta")

#Alinea las secuencias con respecto a virus_seq_not_align
virus_seq_not_align <- OrientNucleotides(virus_seq_not_align)
virus_seq_align <- AlignSeqs(virus_seq_not_align)

#Se guarda el resultado con biostrings
writeXStringSet(virus_seq_align, file="virus_seq_align.fasta")
virus_aligned <- read.alignment("virus_seq_align.fasta", format = "fasta")


#Con seqinr se genera una matriz de distancia donde sombras más oscuras de gris 
#significan una mayor distancia
matriz_distancia <- dist.alignment(virus_aligned, matrix = "similarity")
temp <- as.data.frame(as.matrix(matriz_distancia))

#Creación del árbol con el paquete ape
virus_tree <- nj(matriz_distancia)
class(virus_tree) 
virus_tree <- ladderize(virus_tree)

#Longitud de las secuencias
Lengths <- virus_seq_not_align@ranges@width

#Lee el archivo fasta con bioseq
file <- read_fasta("virus_seqs.fasta")

#Al macena las secuencias en formato DNABin
sequences <- DNA2b(file)

#Concetraciones de nucleotido de las secuencias
Concentraciones <- nuc(sequences)

#Nombres de las secuencias
Nam <- names(virus_sequences)

#Vector con los nombres cuadriplicados
Genoma <- Nams(Nam)

#Etiquetas de los nucletidos
Nucleotidos <- c("A","C","G","T")

#Asigna al eje x las etiquitas de los nucletidos
xlab(Nucleotidos)

#Crea el data frame para ggplot
d <- data.frame(Nucleotidos, Concentraciones, Genoma)

#Crea la grafica de concentraicones
g <- ggplot(data = d, aes(x = Nucleotidos, y = Concentraciones, fill = Genoma))

#Viluación de la grafica de concentraciones
g + geom_bar(stat="identity", position="dodge")

#visualización del arbol
ggtree(virus_tree ) + geom_tiplab()

#Cita 
citation("dplyr")