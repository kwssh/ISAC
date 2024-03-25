function [idx, folder_name] = init_ISAC()

    idx = input('1.save  2.test : ');

    if idx == 1
        folder_name = input('Write foler name : ', 's');
    else
        folder_name = 'qwer';
    end

end