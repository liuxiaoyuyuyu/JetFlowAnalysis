import os
import subprocess

# Directory containing the lists
list_dir = "Run2_tree_list/list_25/"  # Update this to the correct path
name = "data_run2_vn500.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = data_vn500_new.sh
Log        = condor_log/submit_v0.$(Process).log
Output     = condor_log/submit_v0.$(Process).out
Error      = condor_log/submit_v0.$(Process).err
+MaxRuntime = 20000
'''

# Loop over all files in the directory
for filename in os.listdir(list_dir):
    file_path = os.path.join(list_dir, filename)
    if os.path.isfile(file_path):  # Ensure it's a file
        # Add the file path as an argument
        temp = f'''
Arguments  = {file_path}
Queue
        '''
        command_lines += temp

# Write to the JDL file
f.write(command_lines)
f.close()

# Submit the job
subprocess.call(["condor_submit", name])

