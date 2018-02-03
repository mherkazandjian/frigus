clean:
	rm -fvr *.pyc htmlcov

test:
	cd tests && PYTHONPATH=.:${PYTHONPATH} py.test -v --cov --cov-report=html .
