function [varargout]=sr2_find_seq_sets(varargin)
%% [varargout]=sr2_find_seq_sets(varargin)
%
% input: sequence set size S
% output: cell struct with all possible sequence sets of size S matched for
% first-order transitions (first column: training set, second column: test set)
%
% example call:
%               [varargout]=sr2_find_seq_sets(2);
%               [varargout]=sr2_find_seq_sets([2,1]);
%
% --
% gariani@uwo.ca - 2017.10.09

%% global vars
if size(varargin{1},2)>1
    setSize=varargin{1}(1,1);
    triple_set=varargin{1}(1,2);
else
    setSize=varargin{1};
    triple_set=0;
end

% find all possible sequences with each finger only once
allseq=perms([1,2,3,4,5]);

%% remove "runs" (e.g. 1,2,3, ... 5,4,3, ...)
allseq_noRun=zeros(1,size(allseq,2));
rowCount=1;
for r=1:size(allseq,1)
    isRun=zeros(1,size(allseq,2)-2);
    for c=1:size(allseq,2)-2
        if allseq(r,c)==(allseq(r,c+1)-1) && allseq(r,c)==(allseq(r,c+2)-2)
            isRun(1,c)=1;
        elseif allseq(r,c)==(allseq(r,c+1)+1) && allseq(r,c)==(allseq(r,c+2)+2)
            isRun(1,c)=1;
        else
            isRun(1,c)=0;
        end
    end
    if all(isRun==0)
        allseq_noRun(rowCount,:)=allseq(r,:);
        rowCount=rowCount+1;
    else
        continue
    end
end

%% compute all possible combinations of all_seq_noRun
all_comb=nchoosek(1:size(allseq_noRun,1),setSize);
nattempts=size(all_comb,1);
match=cell(1,3); v=1;
for i=1:nattempts
    
    % keep track of n. attempts
    fprintf(1,'Attempt n.%d\n',i);
    
    % pick train set
    all_idx=(1:size(allseq_noRun,1))';
    train_idx=false(size(allseq_noRun,1),1);
    for j=1:size(all_comb,2)
        train_idx=train_idx|(all_idx==all_comb(i,j));
    end
    train_set=allseq_noRun(train_idx,:);
    
    % check first-order transitions in train set
    trans_mat_train=zeros(5,5);
    for r=1:size(train_set,1)
        for c=1:size(train_set,2)-1
            trans_mat_train(train_set(r,c),train_set(r,c+1))=trans_mat_train(train_set(r,c),train_set(r,c+1))+1;
        end
    end
    
    %% find remaining sequences
    remaining_idx=~train_idx;
    remaining_seq=allseq_noRun(remaining_idx,:);
    
    % compute all possible combinations of remaining_seq
    all_comb_remaining=nchoosek(1:size(remaining_seq,1),setSize);
    nattempts_remaining=size(all_comb_remaining,1);
    
    % permute selection until first-order transitions between train and test set match
    for k=1:nattempts_remaining
        
        % pick test set
        test_set=remaining_seq(all_comb_remaining(k,:),:);
        
        % check first-order transitions in test set
        trans_mat_test=zeros(5,5);
        for r=1:size(test_set,1)
            for c=1:size(test_set,2)-1
                trans_mat_test(test_set(r,c),test_set(r,c+1))=trans_mat_test(test_set(r,c),test_set(r,c+1))+1;
            end
        end
        
        % save matching sets
        if isequal(trans_mat_train,trans_mat_test)
            
            if triple_set==1
                
                % find remaining sequences
                remaining_idx_3rd=1:size(remaining_seq,1);
                test_idx=false(size(remaining_idx_3rd,1),1);
                for m=1:size(all_comb_remaining,2)
                    test_idx=test_idx|(remaining_idx_3rd==all_comb_remaining(k,m));
                end
                remaining_seq_3rd=remaining_seq(~test_idx,:);
                
                % compute all possible combinations of remaining_seq_3rd
                all_comb_remaining_3rd=nchoosek(1:size(remaining_seq_3rd,1),setSize);
                nattempts_remaining_3rd=size(all_comb_remaining_3rd,1);
                
                % permute selection until first-order transitions between train and test set match
                for l=1:nattempts_remaining_3rd
                    
                    % pick test set
                    test_set_3rd=remaining_seq_3rd(all_comb_remaining_3rd(l,:),:);
                    
                    % check first-order transitions in test set
                    trans_mat_test_3rd=zeros(5,5);
                    for r=1:size(test_set_3rd,1)
                        for c=1:size(test_set_3rd,2)-1
                            trans_mat_test_3rd(test_set_3rd(r,c),test_set_3rd(r,c+1))=trans_mat_test_3rd(test_set_3rd(r,c),test_set_3rd(r,c+1))+1;
                        end
                    end
                    
                    % save matching sets
                    if isequal(trans_mat_test,trans_mat_test_3rd)
                        
                        fprintf(1,'Match n.%d found!\n\n',v);
                        match{v,1}=train_set;
                        match{v,2}=test_set;
                        match{v,3}=test_set_3rd;
                        v=v+1;
                        
                    else
                        continue
                    end
                end
                
            else
                
                fprintf(1,'Match n.%d found!\n\n',v);
                match{v,1}=train_set;
                match{v,2}=test_set;
                v=v+1;
                
            end
        else
            continue
        end
    end
end

%% save output
%fn='matching_seq_sets.mat';
%save(fn,'match');
varargout{1}=match;

