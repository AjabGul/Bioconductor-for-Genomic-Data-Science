# numerics ----
x <- 1:10
x
# name the each value
names(x) <- letters[1:10]
x
class(x)
# index
x[1:3]
# call multiple columns using concatinate c in single vector
x[c("a", "b")]



# names ----
x <- 1:3
x
names(x) <- c("A", "A", "b")
x
x["A"]



# names_Duplicates------
anyDuplicated(names(x))
names(x) <- c("A", "B", "C")
x
anyDuplicated(names(x))



# numInt------
x <- 1
x
class(x)
x <- 1:3
x
class(x)



# matrices ----

#  creating matrics
x <- matrix(1:9, ncol = 3, nrow = 3)
x
# naming rows
row.names(x) <- c("a", "b", "c")
x
dim(x)
nrow(x)
ncol(x)

# indexing of metrix x[row,col]
x[1:2]
x[1:2,]
x[ ,2:3]
x[1:2,1:2]
x["a",] # first row
x["a",, drop = FALSE]
x[x>5] # subset of matrix


#  create matrix by colmn (by defult)
x <- matrix(1:9, nrow = 3, ncol = 3)
x
#  create matrix by row
x <- matrix(1:9, ncol = 3, nrow = 3, byrow = TRUE)
x






# List ----
# A list can contain mixed data types
my_list <- list(
  name = "Ajab",
  age = 36,
  skills = c("Genomics", "AI", "Cancer Research"),
  married = TRUE
)
my_list
# Access elements
my_list$name       # "Ajab"
my_list$skills[2]  # "AI"
my_list[3]
my_list$skills[3]


# lapply and sapply function on list
x <- list(1:9)
x
lapply(x, mean)        
sapply(x, mean)





# lapplay and sapply ----

# lapply()
# Full form: list apply
# What it does: Applies a function to each element of a list.
# Return type: Always returns a list, no matter what.


# sapply()
# Full form: simplified apply
# What it does: Same as lapply(), but tries to simplify the result.
# Return type: Vector, matrix, or array if possible; otherwise falls back to a list.




#  DataFrames ----

# create a data fram
df <- data.frame(sex = c("M", "M", "F"), age = c(23, 30, 35))
df
df$sex[3]
sapply(df, class)

# convert a data frame into a matrix
as.matrix(df)
# now df/matrix into lists
as.list(df)

# for a complicated type of data we use as from library(methods)
library(methods)
as(df, "matrix")




