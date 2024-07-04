from os import listdir

output = open("results/results.csv", "w")
output.write("correctD,locationD,typeD,wrongD,correctS,locationS,typeS,wrongS,correctR,locationR,typeR,wrongR,correctParalogy,Drate,Lrate,Rrate\n")

for filename in listdir("results"):
    if filename[:8] == "results-" and filename[-4:] == ".txt":
        params = filename[:-4].split("-")
        Drate = float(params[1])
        Lrate = float(params[2])
        Rrate = float(params[3])

        f = open(f'results/{filename}')
        next(f)

        for line in f:
            output.write(f'{line.rstrip()},{Drate},{Lrate},{Rrate}\n')

