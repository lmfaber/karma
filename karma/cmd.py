import subprocess

from logs import logger


class Cmd:
    """ Runs shell commands including pipelines and returns stdout, stderr and the status code. """

    def __init__(self, command):
        self.pipeline = True if "|" in command else False

        if self.pipeline:
            self.cmd = [cmd.strip(" ").split(" ") for cmd in command.split("|")]
        else:
            self.cmd = command.split(" ")

        self.stdout = None
        self.stderr = None
        self.status = None

    def run(self, in_stream=None):
        """ Executes the command(s).
        
        Keyword Arguments:
            in_stream {str} -- Allows to pass a binary string as stdin for the first command. (default: {None})
        """
        PIPE = subprocess.PIPE
        if self.pipeline:
            first_process = self.cmd[0]
            remaining_processes = self.cmd[1:]

            first_process = subprocess.Popen(
                first_process, stdin=PIPE, stdout=PIPE, stderr=PIPE
            )
            stdout, stderr = first_process.communicate(input=in_stream)
            for process in remaining_processes:
                proc = subprocess.Popen(process, stdin=PIPE, stdout=PIPE, stderr=PIPE)
                stdout, stderr = proc.communicate(input=stdout)
        else:
            proc = subprocess.Popen(self.cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            stdout, stderr = proc.communicate(input=in_stream)

        self.stdout = stdout.decode()
        self.stderr = stderr.decode()
        self.status = proc.returncode
        self.__check_command()

    def __check_command(self):
        if self.status != 0:
            logger.error(f"Command failed: {self.cmd}")
            logger.error(self.stderr)
            exit(1)
