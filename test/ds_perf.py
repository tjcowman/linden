import time
import sys
import subprocess
from subprocess import PIPE, run

def main(args):
    test_path = args[1]
    loci_name =  test_path+args[2]
    controls_name =  test_path+args[3]
    cases_name =  test_path+args[4]


    maxThreads = args[5]
    maxUnknown = args[6]

    args = ["--loci", loci_name, "--cases", cases_name, "--controls", controls_name]

    test_bin("./build/linden", args, ["--maxThreads", maxThreads, "--maxUnknown", maxUnknown, "--minMAF", ".05", "--maxMS", "6" ])
    #test_bin("./build/lindenV0", args, ["--maxThreads", maxThreads, "--maxUnknown", maxUnknown, "--minMAF", ".05", "--maxMS", "6"])


def test_bin(bin_path, input_args, other_args):
    t_begin = time.perf_counter()
    result = subprocess.run( [bin_path] + input_args + other_args, stdout=PIPE, text=True ) # stderr=PIPE
    t_end = time.perf_counter()

    data = result.stdout;

    #convert the resulting linden data into a convenient list of lists [id0 = sig, id7 = recip or cutoff]
    tbl = []
    for line in data.split('\n'):
        tbl.append(line.split('\t'))

    total_recip = 0.0
    count_recip = 0;
    total_cutoff = 0.0
    count_cutoff = 0;
    for row in tbl[:-1]: #last line is empty
        if row[7] == 'recip':
            total_recip += float(row[0])
            count_recip += 1
        else:
            total_cutoff += float(row[0])
            count_cutoff += 1

    #print(result.stderr)
    print(t_end-t_begin, total_recip/count_recip, total_cutoff/count_cutoff)


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv)