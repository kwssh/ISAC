function precoder_total = get_precoder_total_letter(num_episode, precoder)

    precoder_total = 0;

    for k = 1:num_episode

        precoder_total = precoder_total + precoder(:,:,k);
    end

end