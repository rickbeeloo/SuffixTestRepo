

with open("sub_test.txt","r") as in_file:
    s = set()
    for line in in_file:
        genome_id = line.split(":")[2].split("_")[0]
        s.add(genome_id)
    print("Elements: ", len(s))
