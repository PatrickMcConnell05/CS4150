import pandas as pd
import math
import statistics

path = "../data/GSE64881_segmentation_at_30000bp.passqc.multibam.txt"

df = pd.read_csv(path, sep='\t')


window_per_np = 2

np_names = df.columns.values[3:]


num_geo_windows = len(df.index) #gets the length of the array of row incicies
num_nps = df.shape[1] - 3 #df.shape returns a tuple of (num_rows, num_cols), we access the cols, minus three



# ------- Evaluate data quality ---------- #

#1.) count # of windows per np, push results to a list

count = 0

arr_windows_per_np = []

for col in df.columns[3:]: #goes through all of the columns
    values = df[col] #assigns a col to values
    count += (values > 0).sum() #sums the values
    arr_windows_per_np.append(count) #adds the 
    count = 0

# print(arr_windows_per_np)

pairing = dict(zip(np_names, arr_windows_per_np)) #key value pairing of the np names and the num of windows present


#2.) get the average windows per np from activity 1

count = 0

for col in df.columns[3:]: #goes through all of the columns
    values = df[col] #assigns a col to values
    count += (values > 0).sum() #sums the values

avg_windows_per_np = count / num_nps

print("Average windows per np:",avg_windows_per_np)


#3.) Calculate the standard deviation
data = pd.Series(arr_windows_per_np)
std_dev = data.std()
print("Standard Deviation:",std_dev)


#4.) Get all nps greater than 2 standard deviations

HIGH_FREQUENCY = avg_windows_per_np + (2 * std_dev) #defined as 2 standard deviations above the average

high_freq_nps = {}
for name, value in pairing.items():
    if value > HIGH_FREQUENCY:
        high_freq_nps[name]=value

# print(len(high_freq_nps))

        


# ------- Estimate radial position of each NP ---------- #


#5.) Take the list of detection frequencies and sort

# arr_windows_per_np.sort()

sorted_pairings = dict(sorted(pairing.items(), key=lambda item: item[1]))
# print("sorted:",sorted_high_freq_nps)


#6.) Split list into five categories

arr_size = len(sorted_pairings)



print("Arr_size:",arr_size)



strongly_apical = []
somewhat_apical = []
neither = []
somewhat_equatorial = []
strong_equatorial = []

radial_group_size = num_nps / 5

pairs = list(sorted_pairings.items())

for index, ((np, frequency)) in enumerate(pairs):
    if index < radial_group_size:
        strongly_apical.append(((np, frequency)))
    elif index < (radial_group_size * 2):
        somewhat_apical.append(np)
    elif index < (radial_group_size * 3):
        neither.append(np)
    elif index < (radial_group_size * 4):
        somewhat_equatorial.append(np)
    elif index < (radial_group_size * 5):
        strong_equatorial.append(np)



# ------- Estimate compaction of each genomic window: ---------- #


comp_group_size = num_nps / 10

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

for index, ((np, frequency)) in enumerate(pairs):
    if index < (comp_group_size):
        one.append(((np, frequency)))
    elif index < (comp_group_size * 2):
        two.append((np, frequency))
    elif index < (comp_group_size * 3):
        three.append((np, frequency))
    elif index < (comp_group_size * 4):
        four.append((np, frequency))
    elif index < (comp_group_size * 5):
        five.append((np, frequency))
    elif index < (comp_group_size * 6):
        six.append((np, frequency))
    elif index < (comp_group_size * 7):
        seven.append((np, frequency))
    elif index < (comp_group_size * 8):
        eight.append((np, frequency))
    elif index < (comp_group_size * 9):
        nine.append((np, frequency))
    elif index < (comp_group_size * 10):
        ten.append((np, frequency))
    
# all_ten = one + two + three + four + five + six + seven + eight + nine + ten
# print("Length of all ten:",len(all_ten))

