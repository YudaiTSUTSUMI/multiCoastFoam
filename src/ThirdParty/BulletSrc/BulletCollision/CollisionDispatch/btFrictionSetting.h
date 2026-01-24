#ifndef FRICTION_SETTING_H
#define FRICTION_SETTING_H

#include <string>
#include <iostream>

enum FrictionCombineMode {
    COMBINE_MULTIPLY,
    COMBINE_MINIMUM,
    COMBINE_AVERAGE
};

class FrictionSetting {
public:
    inline static FrictionCombineMode Mode = COMBINE_MULTIPLY;

    static void SetModeFromString(const std::string& modeStr) {
        if (modeStr == "minimum") {
            Mode = COMBINE_MINIMUM;
            std::cout << "[FrictionSetting] Mode set to MINIMUM\n";
        }
        else if (modeStr == "average") {
            Mode = COMBINE_AVERAGE;
            std::cout << "[FrictionSetting] Mode set to AVERAGE\n";
        }
        else {
            Mode = COMBINE_MULTIPLY;
            std::cout << "[FrictionSetting] Mode set to MULTIPLY (default)\n";
        }
    }
};

#endif
