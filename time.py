import subprocess
import time
import resource

a = [12, 13, 14, 15, 23, 24, 25, 34, 35, 45]

d = ['./sample-data/facebook_combined.snap']

for dd in d:
    for b in a:
        
        with open('test.txt', 'w') as f:
            f.write(b)
        
        with open('test_result_'+str(b)+'.txt', 'w') as f:

            process = subprocess.Popen("./main " + dd, stdout=f)

            process.wait()

            usage = resource.getrusage(resource.RUSAGE_CHILDREN)
            execution_time = usage.ru_utime + usage.ru_stime
            print("Execution Time: %s seconds" % execution_time, stdout=f)

            max_memory = usage.ru_maxrss * resource.getpagesize()
            print("Max Memory Usage: %s bytes" % max_memory, stdout=f)