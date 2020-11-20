#!/usr/bin/python3

def wordfile_to_list(filename):
    """convert word file to a list of word"""
    handler = open(filename, 'r')
    lines = handler.readlines()
    word_list = []
    for line in lines:
        word_list.append(line.rstrip('\n'))
    handler.close()
    return word_list

def wordfile_differences_linearsearch(file1, file2):
    """for each word in the first list, look through the second list for a match"""
    list1 = wordfile_to_list(file1)
    list2 = wordfile_to_list(file2)
    only = []
    for word1 in list1:
        if (not(word1 in list2)):
            only.append(word1)
    return only

def binary_search(sorted_list, element):
    """Search for element in list using binary search.
       Assumes sorted list"""
    # Current active list runs from index_start up to
    # but not including index_end
    index_start = 0
    index_end = len(sorted_list)
    while (index_end - index_start) > 0:
        index_current = (index_end-index_start)//2 + index_start
        if element == sorted_list[index_current]:
            return True
        elif element < sorted_list[index_current]:
            index_end = index_current
        elif element > sorted_list[index_current]:
            index_start = index_current+1
    return False

def wordfile_differences_binarysearch(file1, file2):
    """use binary search to find words that only exist in first file"""
    list1 = wordfile_to_list(file1)
    list2 = wordfile_to_list(file2)
    only = []
    for word1 in list1:
        if(not binary_search(list2,word1)):
            only.append(word1)
    return only

def wordfile_to_dict(filename):
    """read the file and return a list of words"""
    handler = open(filename, 'r')
    lines = handler.readlines()
    word_dict = {}
    for line in lines:
        word_dict[line.rstrip('\n')] = None
    handler.close()
    return word_dict

def wordfile_differences_dictsearch(file1, file2):
    """use dictionary search to find words that only exist in the first file"""
    list1 = wordfile_to_list(file1)
    dict2 = wordfile_to_dict(file2)
    only = []
    for word1 in list1:
        if(not(word1 in dict2)):
            only.append(word1)
    return only