import pandas as pd
import math
import statistics

path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

df = pd.read_csv(path, sep='\t')
df2 = pd.read_csv(path, sep='\t')


# range_to_keep = [21700000, 24100000]

#Creates a deep copy in a dataframe for chrom 13 between the desired coordinates
df2 = df2[(df2['chrom'] == 'chr13') & (df2['start'] >= 21690000) & (df2['stop'] <= 24120000)].copy()


#goes through and eliminates the columns without any windows present
for col in df2.columns[3:]:
    if df2[col].sum() == 0:
        df2 = df2.drop(columns=col)



# print(df2.shape[1] - 3) #minus three to account for the non np cols (num nps)

# ----------------- Basic Statistics -------------------#
#1.) num of genomic windows
num_windows = len(df2.index)
print("Number of genomic windows:",num_windows)

#2.) num nps
num_nps = df2.shape[1] - 3
print("Number of Nuclear Profiles:",num_nps)
# print(df2.head())
# print(df2.tail())

#3.) On average, how many windows are present in an NP?
count = 0

for col in df2.columns[3:]: #goes through all of the columns
    values = df2[col] #assigns a col to values
    count += (values > 0).sum() #sums the values

average = count / num_nps
print("Average Number of Windows in an NP:",average) #printing average number of windows per np

#4.) What is the smallest number of windows present in any NP? The largest?
smallest = 100000000000
largest = 0

for col in df2.columns[3:]:
    values = df2[col] #sets values to be a single column
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


for index, row in df2.iterrows():
    for col in df2.columns[3:]: #goes through every column (element) in the row
        if (row[col] > 0):
            count += 1
            total_count += 1
    if(count < smallest):
        smallest = count
    if(count > largest):
        largest = count
    count = 0


  

average = total_count / (num_windows)
print("Smallest Number of Nps Where a Window is Detected:", smallest)
print("Largest Number of Nps Where a Window is Detected:", largest)
print("Average Number of NPs Where a Window is Detected:", average)

#-------------Find Ranges---------------# 
smallest = 100000000000
largest = 0

for col in df.columns[3:]:
    values = df[col] #sets values to be a single column
    count = (values > 0).sum() #gets the number of 1s in a col
    if(count < smallest):
        smallest = count
    if(count > largest):
        largest = count

# print("Smallest Number of Windows:", smallest)
# print("Largest Number of Windows:", largest)

group_size = (largest - smallest) / 5
# print(group_size)


#6.) work from act2.py --meant to pair names and values up

np_names = df2.columns.values[3:]

df3 = df.loc[:, df.columns.intersection(np_names)].copy() #this creates the dataframe (df3) with the nps that had at least 1 window in the hist1 section
# print(df3)

count = 0

arr_windows_per_np = []


for col in df3.columns: #goes through all of the columns
    values = df3[col] #assigns a col to values
    count += (values > 0).sum() #sums the values
    arr_windows_per_np.append(count) #adds the 
    count = 0

    
pairing = dict(zip(np_names, arr_windows_per_np)) #key value pairing of the np names and the num of windows present
sorted_pairings = dict(sorted(pairing.items(), key=lambda item: item[1]))


# arr_size = len(sorted_pairings)
# print("Arr_size:",arr_size)
# print("Group size:",group_size)


#7.) Split into 5 categories


strongly_apical = []
somewhat_apical = []
neither = []
somewhat_equatorial = []
strong_equatorial = []

pairs = list(sorted_pairings.items())


for np, frequency in pairs:
    if frequency <= group_size:
        strongly_apical.append(((np, frequency)))
    elif frequency <= (group_size * 2):
        somewhat_apical.append((np,frequency))
    elif frequency <= (group_size * 3):
        neither.append((np,frequency))
    elif frequency <= (group_size * 4):
        somewhat_equatorial.append((np,frequency))
    # elif frequency <= (group_size * 5):
    #     strong_equatorial.append((np,frequency))
    else:
        strong_equatorial.append((np,frequency))
        


print("Strongly apical:    ",(strongly_apical))
print("Somewhat apical:    ",(somewhat_apical))
print("Neither:            ",(neither))
print("Somewhat Equatroial:",(somewhat_equatorial))
print("Strongly Equatorial:",(strong_equatorial))

#8.) Split into 10 categories

ten_group_size = (largest - smallest) / 10



one = []
two = []
three = []
four = []
five = []
six = []
seven = []
eight = []
nine = []
ten = []

for np, frequency in pairs:
    if frequency <= (ten_group_size):
        one.append(((np, frequency)))
    elif frequency <= (ten_group_size * 2):
        two.append((np, frequency))
    elif frequency <= (ten_group_size * 3):
        three.append((np, frequency))
    elif frequency <= (ten_group_size * 4):
        four.append((np, frequency))
    elif frequency <= (ten_group_size * 5):
        five.append((np, frequency))
    elif frequency <= (ten_group_size * 6):
        six.append((np, frequency))
    elif frequency <= (ten_group_size * 7):
        seven.append((np, frequency))
    elif frequency <= (ten_group_size * 8):
        eight.append((np, frequency))
    elif frequency <= (ten_group_size * 9):
        nine.append((np, frequency))
    # elif frequency <= (ten_group_size * 10):
    #     ten.append((np, frequency))
    else:
        ten.append((np, frequency))

print("Category 1:",len(one))
print("Category 2:",len(two))
print("Category 3:",len(three))
print("Category 4:",len(four))
print("Category 5:",len(five))
print("Category 6:",len(six))
print("Category 7:",len(seven))
print("Category 8:",len(eight))
print("Category 9:",len(nine))
print("Category 10:",len(ten))
