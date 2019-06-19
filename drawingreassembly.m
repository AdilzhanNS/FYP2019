

figure;
subplot(2,2,1);
imshow(p);
hold on;
plot(subsetp_1(:,2), subsetp_1(:,1), 'g*', 'LineWidth',1);
hold off;
subplot(2,2,2);
imshow(q);
hold on;
plot(subsetq_2(:,2), subsetq_2(:,1), 'r*', 'LineWidth',1);
plot(subsetq_3(:,2), subsetq_3(:,1), 'b*', 'LineWidth',1);
hold off;
subplot(2,2,3);
imshow('simple3.jpg');
hold on;
plot(subsetw_1(:,2), subsetw_1(:,1), 'g*', 'LineWidth',1);
plot(subsetw_2(:,2), subsetw_2(:,1), 'r*', 'LineWidth',1);
hold off;
subplot(2,2,4);
imshow(x);
hold on;
plot(subsetx_3(:,2), subsetx_3(:,1), 'b*', 'LineWidth',1);
hold off;

%and so on for all the fragments presented
%example of the way the reassembly is done
figure;
subplot(1,3,1);
imshow('plate_1.png');
hold on;
plot(subsetp_1(:,2),subsetp_1(:,1), 'g*', 'LineWidth',1);
plot(substep_2(:,2),substep_2(:,1), 'r*', 'LineWidth',1);
hold off;
subplot(1,3,2);
imshow('plate_2.png');
hold on;
plot(subsetq_1(:,2),subsetq_1(:,1), 'g*', 'LineWidth',1);
hold off;
subplot(1,3,3);
imshow('plate_3.png');
hold on;
plot(subsetw_2(:,2),subsetw_2(:,1), 'r*', 'LineWidth',1);
hold off;




