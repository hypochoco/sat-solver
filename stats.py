
# process information from logs

def log_to_values(path): # log -> list of numbers
    output = []
    with open(path) as f:
        for line in f:
            literal = line.split(", ")
            time_str = literal[0]
            output.append(float(time_str[6:]))
    return output

if __name__ == "__main__":
    output = log_to_values("logs/3sat_dpll.log")
    print(output)