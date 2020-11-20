#!/usr/bin/python3
import timeit
start_time = timeit.default_timer()

import handin6
differences = handin6.wordfile_differences_binarysearch("british-english", "american-english")

time_spent = timeit.default_timer() - start_time
