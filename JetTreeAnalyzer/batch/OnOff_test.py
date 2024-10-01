import math
import subprocess
name = "data_run3_vn500_test.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = data_vn500_test.sh
Arguments  = 000
Log        = condor_log/submit_v0.$(Process).log
Output     = condor_log/submit_v0.$(Process).out
Error      = condor_log/submit_v0.$(Process).err
+MaxRuntime =3000
Queue
'''

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
