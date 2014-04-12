#include <stdio.h>

#include <gmock/gmock.h>

GTEST_API_ int main(int argc, char **argv) {
	printf("Running NEB acceptance tests\n");

	unsigned int oldflags;
	//_controlfp_s(&oldflags, ~(_SW_OVERFLOW | _SW_ZERODIVIDE | _SW_INVALID), _MCW_EM);

	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
