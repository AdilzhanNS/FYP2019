function call_prompt
prompt = 'Do you want more contours? Y/N [Y]: ';
str = input(prompt,'s');
if isempty(str)
    str = 'Y';
end
if (str=='Y')
    f = msgbox('Run the extract_COI algorithm again');
end

if (str=='N')
    f = msgbox('Proceed with other operations');
end

end
    