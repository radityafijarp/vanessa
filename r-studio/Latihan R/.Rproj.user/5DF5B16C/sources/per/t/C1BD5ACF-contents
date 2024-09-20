# Load the mtcars 
mtcars_data=as.data.frame(mtcars)

##Exercise 1.1
print(head(mtcars_data,n=3))

##Exercise 1.2
#a
filtered_5_gears=subset(mtcars_data, gear == 5)
print(filtered_5_gears)

#b (total rows)
total_car_5=dim(filtered_5_gears)[1]
print(total_car_5)

#c (Extract only the names of the cars that have 5 gears, as a character vector.)
row_filtered_vector <- rownames(filtered_5_gears)
print(row_filtered_vector)

#d Display the subset of the table with only cars that have 4 gears and at least 100 hp.
filtered_5_gears_100_hp=subset(subset(mtcars_data, gear == 5),hp>=100)
print(filtered_5_gears_100_hp)

##1.3
#a
toyota_rows <- mtcars_data[grep("Toyota", rownames(mtcars_data)), ]
print(toyota_rows)

#b
mercedes_rows=mtcars_data[grep("Merc", rownames(mtcars_data)), ]
mercedes_rows= mercedes_rows[, c("mpg", "cyl")]
print(mercedes_rows)


library(ggplot2)