#!/bin/bash
# Copyright (C) Matthias Kesenheimer - All Rights Reserved
# Written by Matthias Kesenheimer <m.kesenheimer@gmx.net>, 2017

make pastegnudata
cd ./chaIchaJ && make -j4 all
cd ../neuIchaJ && make -j4 all
cd ../neuIneuJ && make -j4 all
