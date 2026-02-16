import pandas as pd

path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

df = pd.read_csv(path, sep='\t')


#1.) Number of genomic windows
num_geo_windows = len(df.index) #gets the length of the array of row incicies
print("Number of Genomic Windows:",num_geo_windows) #gets the numer of rows, minus the column row


#2.) Number of NPs
num_nps = df.shape[1] - 3 #df.shape returns a tuple of (num_rows, num_cols), we access the cols, minus three
print("Number of Nuclear Profiles (NP):",num_nps)


#3.) On average, how many windows are present in an NP?
count = 0

for col in df.columns[3:]: #goes through all of the columns
    values = df[col] #assigns a col to values
    count += (values > 0).sum() #sums the values
    

average = count / num_nps
print("Average Number of Windows in an NP:",average) #printing average number of windows per np

#4.) What is the smallest number of windows present in any NP? The largest?
smallest = 100000000000
largest = 0

for col in df.columns[3:]:
    values = df[col] #sets values to be a single column
    count = (values > 0).sum() #gets the number of 1s in a col
    if(count < smallest):
        smallest = count
    if(count > largest):
        largest = count

print("Smallest Number of Windows:", smallest)
print("Largest Number of Windows:", largest)

#5.) On average, what is the number of NPs in which a window is detected? The smallest? The largest?
smallest = 100000000
largest = 0

count = 0 #used for the smallest and largest number of nps where a window is detected
total_count = 0 #for the total number of nps


for index, row in df.iterrows():
    for col in df.columns[3:]: #goes through every column (element) in the row
        if (row[col] > 0):
            count += 1
            total_count += 1
    if(count < smallest):
        smallest = count
    if(count > largest):
        largest = count
    count = 0


  

average = total_count / (num_geo_windows)
print("Smallest:", smallest)
print("Largest:", largest)
print("Average Number of NPs Where a Window is Detected:", average)