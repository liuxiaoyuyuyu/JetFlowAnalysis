import math
import subprocess
name = "data_run3_vn500.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = data_vn500.sh
Arguments  = 000
Log        = condor_log/submit_v0.$(Process).log
Output     = condor_log/submit_v0.$(Process).out
Error      = condor_log/submit_v0.$(Process).err
+MaxRuntime =20000
Queue
'''

with open('MCList_Run3_2023_D.list') as fs:
    line_count = 0
    for line in fs:
        line_count += 1

line_count_frac = line_count/25.0 #for all pthatpythia
file_count = math.ceil(line_count_frac)

for i in range(1, int(file_count)):
   temp = '''
Arguments  = %03d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
