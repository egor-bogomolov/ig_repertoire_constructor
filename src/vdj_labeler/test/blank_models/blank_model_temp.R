# Function to take columns (in @id.filename and @probs.filename)
# from matlab and print model in needed form to @model.filename.

get.probabilities.for.id <- function(id.filename, probs.filename, model.filename) {
  df <- read.csv(id.filename, sep='\t', header=F)
  probs <- read.csv(probs.filename, sep='\t', header=F)
  probs <- as.vector(unlist(probs))
  df$V3 <- format(probs[df$V2], scientific=FALSE)
  df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
  write.table(df[, c(1,3)], file=model.filename, sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  NULL
}

## genV
## Getting probabilities for V-genes.
df <- read.csv("V_id/V_id.csv", sep='\t', header=F)
probs <- read.csv("V_id/V_probs.csv", sep='\t', header=F)
probs <- as.vector(unlist(probs))
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="V_id/V_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

## genD
df <- read.csv("blank_model_temp.csv", sep='\t', header=F)
probs <- read.csv("blank_model_temp2.csv", sep='\t', header=F)
probs <- rowSums(probs)
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="blank_model_temp3.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

## genJ
df <- read.csv("blank_model_temp.csv", sep='\t', header=F)
probs <- read.csv("blank_model_temp2.csv", sep='\t', header=F)
probs <- colSums(probs)
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="blank_model_temp3.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# del
df <- read.csv("blank_model_temp.csv", sep='\t', header=F)
probs <- as.matrix(t(read.csv("blank_model_temp2.csv", sep='\t', header=F)))
probs <- apply(probs, c(1, 2), function(x) format(x, scientific=FALSE, digits = 15))
df <- cbind(df, probs[df$V2,])
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, -2], file="blank_model_temp3.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )


# del D genes
df <- read.csv("blank_model_temp.csv", sep='\t', header=F)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
temp <- list.files(path = "blank_models/Ddels", pattern="*.csv")
myfiles <- lapply(temp, function(x) { read.csv(paste("blank_models/Ddels/", x, sep = ''), header=F, sep=',') } )
df$V3 <- myfiles[df$V2]
df$Left <- lapply(df$V3, rowSums)
df$Right <- lapply(df$V3, colSums)
df$Left <- lapply(df$Left, function(x) format(x, scientific = FALSE))
df$Right <- lapply(df$Right, function(x) format(x, scientific = FALSE))
df$Left <- sapply(df$Left, unlist)
df$Right <- sapply(df$Right, unlist)

write.table(cbind(df[, 1], t(simplify2array(df[, 4]))), file="delDLeft.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
