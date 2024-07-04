import os
import getopt
import sys
import subprocess

spname = "species"
gname = "gene"

#default values
nspecies = 10

replicates = 100

results = open("results.txt", "w")
results.write("correctD,locationD,typeD,wrongD,correctS,locationS,typeS,wrongS,correctR,locationR,typeR,wrongR,correctParalogy\n")


#process options
opts, args = getopt.getopt(sys.argv[1:], "s:d:l:r:",["replicates=","seed="])

simgene_com = ["./simulate_gene", spname]

for opt, arg in opts:
    if opt == "-s":
        nspecies = int(arg)
    elif opt == "-d":
        simgene_com.append("--Drate")
        simgene_com.append(arg)
    elif opt == "-l":
        simgene_com.append("--Lrate")
        simgene_com.append(arg)
    elif opt == "-r":
        simgene_com.append("--Rrate")
        simgene_com.append(arg)
    elif opt == "--replicates":
        replicates = int(arg)
    elif opt == "--seed":   #does not work
        seed = int(arg)

os.system("ulimit -s unlimited")

more = 0

for r in range(replicates+more):
    os.system("./simulate_species " + str(nspecies) + " > " + spname)
    try:
        p = subprocess.run(simgene_com, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError:
        more += 1
        continue

    simgene = p.stdout.split("\n")

    gene_out = open(gname, "w")

    true_recon = []
    true_para = []

    simgene_iter = iter(simgene)

    while (l := next(simgene_iter, "end")) != "end":
        if l.rstrip() == "Gene network":
            while (l := next(simgene_iter).rstrip()):
                gene_out.write(l + "\n")

        if l.rstrip() == "Reconciliation":
            while (l := next(simgene_iter).rstrip()):
                true_recon.append(l.split()[1:])

        if l.rstrip() == "Gene paralogy information":
            while (l := next(simgene_iter).rstrip()):
                true_para.append(l.split())

    gene_out.close()

    mpr_recon = []
    mpr_para = []

    #perecon = os.popen("./perecon " + gname + " " + spname)
    try:
        p = subprocess.run(["./perecon",gname,spname], capture_output=True, text=True, check=True, timeout=600)
    except subprocess.CalledProcessError:
        more += 1
        continue
    except subprocess.TimeoutExpired:
        more += 1
        continue

    perecon = p.stdout.split("\n")

    written = False

    perecon_iter = iter(perecon)

    while (l := next(perecon_iter, "end")) != "end":
        if l.rstrip() == "MPR":
            while (l := next(perecon_iter).rstrip()):
                mpr_recon.append(l.split()[1:])
            written = True

        if l.rstrip() == "Gene paralogy information":
            while (l := next(perecon_iter).rstrip()):
                mpr_para.append(l.split())

    #print(true_recon)
    #print(mpr_recon)

    events = {"D":0,"S":1,"R":2}
    counts = [[0]*4,[0]*4,[0]*4]    #first index: D, S, R; second index: correct location and type, correct location only, correct type only, incorrect

    for i,n in enumerate(true_recon):
        if n[1] in events:
            #print(n, events[n[1]], 0 if n == mpr_recon[i] else 1)
            counts[events[n[1]]][(0 if n[0] == mpr_recon[i][0] else 2) + (0 if n[1] == mpr_recon[i][1] else 1)] += 1

    #process paralogy information
    correct_para = 0;
    total_para = 0;
    for p in true_para:
        for mp in mpr_para:
            if {p[0], p[1]} == {mp[0], mp[1]} and (p[2],p[3]) == (mp[2],mp[3]):
                if p[4] == mp[4]:
                    correct_para += float(p[3]) - float(p[2]);
                total_para += float(p[3]) - float(p[2]);

    if total_para == 0:
        correct_para = 1.0
        total_para = 1.0

    results.write(','.join([','.join(map(str,c)) for c in counts]) + "," + str(correct_para/total_para) + "\n")

#what do i need to do?
#replicate  ./
#change # species   ./
#bug-fixing ./
#change Drate/Lrate/Rrate for genes  ./
#finer-grained comparison: correct placement, correct event type    ./
#compare paralogy   ./
