def joiner(file1, file2, path_to_output="NoNameGiven"):
    """
    takes in two textfiles and joins them for better comparison
    
    input:
    file*: the files you want to compare
    """
    from pathlib import Path
    # taking in data and making it accessable
    data1 = open(file1, "r")
    data2 = open(file2, "r")
    lines1 = data1.readlines()
    lines2 = data2.readlines()
    total = min(len(lines1), len(lines2))

    # getting output
    if filepath == "NoNameGiven":
        raise ValueError("no path specified")
    Path(path_to_output).touch()
    output_file = open(path_to_output, "w")
    header = f"cuhf         uhf\n"
    output_file.writelines(header)
    for line in range(total):
        output_file.writelines(f"{lines1[line]:>12} {lines2[line]:>12}\n")
        print(f"{lines1[line]:>12} {lines2[line]:>12}")
    output_file.close()

    joiner("Bachelor_project/Bachelor_project/datafiles/h20_cuhf.txt", "Bachelor_project/Bachelor_project/datafiles/h20_uhf.txt", path_to_output="Bachelor_project/Bachelor_project/datafiles/compare.txt")

