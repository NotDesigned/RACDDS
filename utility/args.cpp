//
// Created by yy on 11/30/23.
//

#include "args.h"


Args::Args() {
    args_["-path"] = "./";
    args_["-t"] = "u";
    args_["-a"] = "e";
    args_["-eps"] = "0";
    args_["-red"];
    args_["-alloc"];
    args_["-ext"];
    args_["-ver"];
    args_["-o"] = "0";
    args_["-seq"] = "f";
    args_["-vw"] = "f";
    args_["-gamma"] = "1";
//    args_["-p"] = "1";
    args_["-exp"] = "t";
    args_["-it"] = "-1";
    args_["-dc"] = "t";
    args_["-ra"] = "f";
    args_["-stable"] = "f";
    args_["-map"] = "t";
    args_["-res"] = "f";
    args_["-width"] = "2";
    args_["-multi"] = "f";
    args_["-estable"] = "t";
    args_["-stats"] = "t";
    args_["-sample"] = "f";
    args_["-rate"] = "1";
    args_["-printc"] = "f";
    args_["-wshrink"] = "f";
    args_["-initwcore"] = "f";
    args_["-adam"] = "f";
    args_["-original"] = "f";
    args_["-log"] = "f";
}


void Args::argsParse(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (std::find(options_.begin(), options_.end(), argv[i]) == options_.end()) {
            //error
            printf("Parsing Error on %s\n", argv[i]);
            exit(-1);
        }
        if (i + 1 > argc) {
            //error
            exit(-1);
        }
//        printf("%s: %s\n", argv[i], argv[i+1]);
        args_[argv[i]] = argv[i + 1];
        i++;
    }
}

std::string Args::getOption(const std::string &option) {
    if (std::find(options_.begin(), options_.end(), option) == options_.end()) {
        //error
        std::cout <<"GetOption Error on "<< option <<std::endl;
        exit(-1);
    }
    return args_[option];
}
