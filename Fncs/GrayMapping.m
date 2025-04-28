
function const = GrayMapping(M, constType)
% Gray Mapping for digital modulations.
%
% Parameters
% ----------
% M : int
%     modulation order
% constType : 'qam', 'psk', 'pam' or 'ook'.
%     type of constellation.
%
% Returns
% -------
% const : list
%     list with constellation symbols (sorted according their corresponding
%     Gray bit sequence as integer decimal).

if M ~= 2 && constType == "ook"
    warning("OOK has only 2 symbols, but M != 2. Changing M to 2.")
    M = 2;
end


if constType=="pam" || constType=="ook" 
    L = (M - 1);
else 
    L=sqrt(M) - 1;
end

bitsSymb = log2(M);

code = GrayCode(bitsSymb); % use a custom function to generate gray code
if constType == "ook"
    const = 0:1;
elseif constType == "pam"
    const = -L:2:L;
elseif constType == "qam"
    PAM = -L:2:L;
    PAM = reshape(PAM, [1, length(PAM)]);

    % generate complex square M-QAM constellation
    const = repmat(PAM, [L + 1, 1]);
    const = const + 1i * flipud(const.');
    
    for ind = 1:2:L+1
        const(ind, :) = flip(const(ind, :), 2);
    end
elseif constType == "psk"
    pskPhases = 0:2 * pi / M:2 * pi;

    % generate complex M-PSK constellation
    const = exp(1i * pskPhases);
end
const = reshape(const, [M, 1]);
const_ = zeros(M, 2);

for ind = 1:M
    const_(ind, 1) = const(ind, 1); % complex constellation symbol
    const_(ind, 2) = bin2dec(cell2mat(code(ind))); % mapped bit sequence (as integer decimal)
end
% sort complex symbols column according to their mapped bit sequence (as integer decimal)
const = sortrows(const_, 2);
const = const(:, 1);

if constType=="pam" || constType=="ook"
    const = real(const);
end

end
