# Function to print attributes available extracted
# from a VCF file. It takes the 'attributes' output 
# from VCFtoTable function. 

def attribute_guide(attribute_dict):
    for key in attribute_dict:
        print("\n**** %s Attributes ****\n" % key)
        for attribute in attribute_dict[key]:
            print("%s:\nDescription: %s\nType: %s\n" % 
                  (attribute, attribute_dict[key][attribute]['Description'], 
                   attribute_dict[key][attribute]['Type']))
