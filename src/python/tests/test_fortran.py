import py.test
from subprocess import Popen

def run_test(test_name):

    process = Popen(test_name)
    process.wait()
    if process.returncode > 0:
        raise ValueError('test failed')

def test_1():
    py.test.skip('not implemented yet')
    run_test('/path/to/cooling_function/src/forantan/test1.ext')

def test_2():
    py.test.skip('not implemented yet')
    run_test('/path/to/cooling_function/src/forantan/test2.ext')
