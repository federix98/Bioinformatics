if __name__ == "__main__":
    files = ["hairpin", "mature"]
    for file in files:
        print("doing", file)
        with open("%s_homo.fa" % file, "w") as out:
            with open("%s.fa" % file, 'r') as file:
                row = next(file)
                while row:
                    try:
                        if ">" in row and "Homo" in row:
                            out.write(row)
                            print(row)
                            row = next(file)
                            while not ">" in row:
                                print(row)
                                out.write(row)
                                row = next(file)
                        else:
                            row = next(file)
                    except Exception:
                        row = None
                        print("Exception")