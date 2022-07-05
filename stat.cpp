#include "stat.h"

namespace snu {
	
	void initStat(StatResult& result) {
		result.basicstat = false;
		result.connectstat = false;
        result.biconnectedstat = false;
		result.countstat = false;
	}
}
