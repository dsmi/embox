function suite_all
% suite_all : testsuite including all the solver components tests.
%

addTestSuite('suite_common');
addTestSuite('suite_mesh');
addTestSuite('suite_tlines');
addTestSuite('suite_solver');
