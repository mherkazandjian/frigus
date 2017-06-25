clean:
	rm -fvr *.pyc

test:
	cd tests && PYTHONPATH=..:$PYTHONPATH py.test -v --cov --cov-report=html .
