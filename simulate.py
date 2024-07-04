import os

os.system("ulimit -s unlimited")

for rate in range(11):
    os.system(f'python3 wrapper.py -d {rate/10.0} -l {rate/10.0} --replicates 1000')
    os.system(f'mv results.txt new_results/results-{rate/10.0}-{rate/10.0}-0.5.txt')

for rate in range(11):
    os.system(f'python3 wrapper.py -d {rate/10.0} -r {rate/10.0} --replicates 1000')
    os.system(f'mv results.txt new_results/results-{rate/10.0}-0.5-{rate/10.0}.txt')
