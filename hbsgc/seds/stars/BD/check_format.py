from os import listdir

for a_file in listdir():
    with open(a_file) as f:
        for line in f:
            if 'E' not in line:
                print(a_file)
