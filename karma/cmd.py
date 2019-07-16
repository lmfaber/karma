import subprocess
from logs import logger

class Cmd():

    """
    Runs shell commands including pipelines and returns stdout, stderr and the status code.
    """

    def __init__(self, command):
        self.pipeline = True if '|' in command else False

        if self.pipeline:
            self.cmd = [cmd.strip(' ').split(' ') for cmd in command.split('|')]
        else:
            self.cmd = command.split(' ')

        self.stdout = None
        self.stderr = None
        self.status = None
    
    def run(self):
        """
        Runs the command(s).
        """
        PIPE = subprocess.PIPE
        if self.pipeline:
            first_process = self.cmd[0]
            remaining_processes = self.cmd[1:]

            proc_1 = subprocess.Popen(first_process, stdout=PIPE, stderr=PIPE)
            for process in remaining_processes:
                proc = subprocess.Popen(process, stdin=proc_1.stdout, stdout=PIPE, stderr=PIPE)
                proc_1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
                proc_1 = proc     
        else:
            proc = subprocess.Popen(self.cmd, stdout=PIPE, stderr=PIPE)

        stdout, stderr = proc.communicate()
        self.stdout = stdout.decode()
        self.stderr = stderr.decode()
        self.status = proc.returncode
