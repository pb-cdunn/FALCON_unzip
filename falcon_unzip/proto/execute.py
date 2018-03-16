#! /usr/bin/env python2.7

# import glob
import traceback
import subprocess
from time import gmtime, strftime

"""
Logs a message to a given file handler. Adds the current date and
time to the front.
"""
def log(message, fp_log, silent=False):
  if (fp_log == None):
    return
  if (silent == True):
    return

  timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
  if (message != ''):
    prefix = '[%s] ' % (timestamp)
  else:
    prefix = ''

  fp_log.write('%s%s\n' % (prefix, message))
  fp_log.flush()

"""
Executes a given shell command and logs the attempt.
"""
def execute_command(command, fp_log, dry_run=True):
  if (dry_run == True):
      log('Executing (dryrun): "%s".' % (command), fp_log);
      log('\n', fp_log)
      return 0
  if (dry_run == False):
      log('Executing: "%s".' % (command), fp_log);
      rc = subprocess.call(command, shell=True)
      if (rc != 0):
          log('ERROR: subprocess call returned error code: %d.' % (rc), fp_log);
          log('Traceback:', fp_log);
          traceback.print_stack(fp_log)
          exit(1)
      return rc


"""
Executes a given shell command and logs the attempt, but also gets
the stdout, stderr and the returne code of the program.
"""
# def execute_command_with_ret(dry_run, command, silent=False):
def execute_command_with_ret(command, fp_log, dry_run=True):
  # if (silent == False):
  #   sys.stderr.write('Executing command: "%s"\n' % command);
  log('Executing (dryrun): "%s".' % (command), fp_log);
  log('\n', fp_log)
  if (dry_run == False):
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    [output, err] = p.communicate()
    rc = p.returncode

  log('\n', fp_log)

  return [rc, output, err]
