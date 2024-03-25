import math
import subprocess
name = "approv_data_vn260.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = approv_data/data_vn260.sh
Arguments  = 000
Log        = approv_data/err_log_out260/submit_v0.$(Process).log
Output     = approv_data/err_log_out260/submit_v0.$(Process).out
Error      = approv_data/err_log_out260/submit_v0.$(Process).err
+MaxRuntime =3000
Queue
'''

with open('raw_mini.txt') as fs:
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
