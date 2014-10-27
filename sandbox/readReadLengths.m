function [T_x, T_hist] = readReadLengths()

filename = '../raw_data/HL7_insert_size_plain.dat';
readOut = dlmread(filename, '\t', 1,0);
T_x  = readOut(:,1);
T_hist  = readOut(:,2);