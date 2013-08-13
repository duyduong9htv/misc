function out = string_pad(string_in, desired_length)
%function out = string_pad(string_in, desired_length)
% pads to the front of a stirng with zeros until length of string is desired. 
%: ex: string_pad('5', 6)
	out = string_in; 
    while length(out) < desired_length
        out = ['0' out];
	end 
end
