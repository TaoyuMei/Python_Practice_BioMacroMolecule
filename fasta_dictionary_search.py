#!/usr/bin/python3

def read_fasta(file):
    """read fasta file and create dictionary"""
    import re
    handler=open(file,'r')
    lines=handler.readlines()
    dictionary={}
    value=""
    key=""
    for line in lines:
        if line[0]==">":
            if(value!="" and key!=""):
                dictionary[key]=value
            pattern=re.compile('>(.+)\n')
            key=pattern.match(line).group(1)
            value=""
        else:
            value+=line.rstrip('\n')
    dictionary[key]=value
    handler.close()
    return dictionary

def find_prot(dictionary,prot_name):
    """search for a given protein and find its sequence"""
    if dictionary.get(prot_name)!=None:
        return dictionary[prot_name]
    else:
        print("Error: protein name "+prot_name+" not found")
        return None

def find_prot2(dictionary,re_expr):
    """return all the keys in the dictionary that match the regular expression"""
    import re
    list_keys=[]
    pattern=re.compile(re_expr)
    keys=dictionary.keys()
    for key in keys:
        match=pattern.search(key)
        if match:
            list_keys.append(match.group())
    return list_keys