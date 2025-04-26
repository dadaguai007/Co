function [u1x,u1y] = DDLMS_v0(u0x,u0y,mu,ntaps,sps,P)
    L = floor((length(u0x)-ntaps)/sps);
    Lp = floor(L/2);
    err = nan(L,2);
    h = zeros(ntaps,2,2);
    mid = ceil((ntaps+1)/2);
    h_init = 1;
    h(mid,:,:) = diag([h_init h_init]);
    hxx = squeeze(h(:,1,1));
    hxy = squeeze(h(:,1,2));
    hyx = squeeze(h(:,2,1));
    hyy = squeeze(h(:,2,2));
    u_hat = nan(L,2);
    u_temp = nan(ntaps,2);
    temp = 1:ntaps;
    field_tp = [u0x, u0y];
    for i = 1:Lp
        idx = temp + (i-1)*sps;
        field = field_tp(idx,:);
        u_temp(:,1) = sum([hxx, hxy].*field, 2); % [ntaps,1]
        u_temp(:,2) = sum([hyx, hyy].*field, 2); % [ntaps,1]
        dn = sum(u_temp,1); % [1,2]
        mu_mat = repmat(mu,[2,2]);
        [dec_x, ~] = Decision(dn(1),P);
        [dec_y, ~] = Decision(dn(2),P);
        % 误差计算
        err(i,1) = dec_x - dn(1);
        err(i,2) = dec_y - dn(2);
        % 抽头更新
        hxx = hxx + mu_mat(1) * err(i,1) * conj(field_tp(idx,1));
        hyx = hyx + mu_mat(2) * err(i,2) * conj(field_tp(idx,1));
        hxy = hxy + mu_mat(3) * err(i,1) * conj(field_tp(idx,2));
        hyy = hyy + mu_mat(4) * err(i,2) * conj(field_tp(idx,2));
    end

    for k = 1:L
        idx = temp + (k-1)*sps;
        field = field_tp(idx,:);
        u_temp(:,1) = sum([hxx, hxy].*field, 2);
        u_temp(:,2) = sum([hyx, hyy].*field, 2);
        u_hat(k,:) = sum(u_temp,1);
        dn = u_hat(k,:);
        mu_mat = repmat(mu,[2,2]);
        [dec_x, ~] = Decision(dn(1),P);
        [dec_y, ~] = Decision(dn(2),P);
        % 误差计算
        err(k,1) = dec_x- dn(1);
        err(k,2) = dec_y- dn(2);
        % 抽头更新
        hxx = hxx + mu_mat(1) * err(k,1) * conj(field_tp(idx,1));
        hyx = hyx + mu_mat(2) * err(k,2) * conj(field_tp(idx,1));
        hxy = hxy + mu_mat(3) * err(k,1) * conj(field_tp(idx,2));
        hyy = hyy + mu_mat(4) * err(k,2) * conj(field_tp(idx,2));

    end
    u1x = u_hat(:,1);
    u1y = u_hat(:,2);
  
end
function [dec,index] = Decision(u, P)
    u = repmat(u,[length(P),1]);
    [~,index] = min(abs(u-P));
    dec = P(index);
end