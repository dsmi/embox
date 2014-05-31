function test_wrapidx
%

widx = wrapidx([ -2 -1 0 1 5 ], 4);

assertEquals([ 2 3 4 1 1 ], widx);
