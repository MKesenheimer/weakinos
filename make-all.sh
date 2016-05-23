#!/bin/bash

make pastegnudata
cd ./chaIchaJ && make -j4 all
cd ../neuIchaJ && make -j4 all
cd ../neuIneuJ && make -j4 all
