function plot_residuals(constellations, RES_PHASE, RES_CODE, outliers_PHASE, outliers_CODE, filerootOUT)

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2017 Mirko Reguzzoni, Eugenio Realini
%  Written by:
%  Contributors:     ...
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------


if length(find(isfinite(RES_PHASE))) > 1
    plot_phase=1;
else
    plot_phase=0;
end

if length(find(isfinite(RES_CODE))) > 1
    plot_code=1;
else
    plot_code=0;
end

outliers_PHASE(isnan(outliers_PHASE)) = 0;
outliers_CODE(isnan(outliers_CODE)) = 0;

nrows = 7;
ncols = 6;

systems = fieldnames(constellations);
for s = 1 : 6 %numel(unique(constellations.systems))
    if constellations.(systems{s}).enabled == 1

        %PHASE GRAPHS (no outliers)
        if plot_phase == 1
            f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
            paperSize = get(f,'PaperSize');
            set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
            for i=1:constellations.(systems{s}).numSat

                figure(f);
                subplot(nrows,ncols,i);
                hold on;
                title(sprintf('%c%02d',constellations.(systems{s}).sysID,constellations.(systems{s}).PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
                grid on;
                set(gca,'FontName','Verdana');
                set(gca,'FontSize',7);
                ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                plot(RES_PHASE(constellations.(systems{s}).indexes(i),~outliers_PHASE(constellations.(systems{s}).indexes(i),:))*1000,'.b');
           end

            try
                print(f , '-dpdf', [filerootOUT '_' systems{s} '_PHASE_residuals']);
            catch ex
                log = Logger.getInstance();
                log.addError(sprintf('saving phase residuals PDF failed - %s\n', ex.message))
            end
            %remove figure
            close(f);
        end

        %CODE GRAPHS (no outliers)
        if plot_code == 1
            f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
            paperSize = get(f1,'PaperSize');
            set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);

            for i=1:constellations.(systems{s}).numSat

                figure(f1);
                subplot(nrows,ncols,i);
                hold on;
                title(sprintf('%c%02d',constellations.(systems{s}).sysID,constellations.(systems{s}).PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
                grid on;
                set(gca,'FontName','Verdana');
                set(gca,'FontSize',7);
                ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                plot(RES_CODE(constellations.(systems{s}).indexes(i),~outliers_CODE(constellations.(systems{s}).indexes(i),:)),'.b');
            end
            try
                print(f1 , '-dpdf', [filerootOUT '_' systems{s} '_CODE_residuals']);
            catch ex
                log = Logger.getInstance();
                log.addError(sprintf('saving code residuals PDF failed - %s\n', ex.message))
            end


            close(f1);
        end

        %PHASE GRAPHS (with outliers)
        if plot_phase == 1
            f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
            paperSize = get(f,'PaperSize');
            set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
            for i=1:constellations.(systems{s}).numSat

                figure(f);
                subplot(nrows,ncols,i);
                hold on;
                title(sprintf('%c%02d',constellations.(systems{s}).sysID,constellations.(systems{s}).PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
                grid on;
                set(gca,'FontName','Verdana');
                set(gca,'FontSize',7);
                ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                plot(RES_PHASE(constellations.(systems{s}).indexes(i),:)*1000,'.b');
                hold on
                plot(find(outliers_PHASE(constellations.(systems{s}).indexes(i),:) == 1),RES_PHASE(constellations.(systems{s}).indexes(i), outliers_PHASE(constellations.(systems{s}).indexes(i),:) == 1)*1000,'.r');
            end
            try
                print(f , '-dpdf', [filerootOUT '_' systems{s} '_PHASE_residuals_with_OUTLIERS']);
            catch ex
                log = Logger.getInstance();
                log.addError(sprintf('saving phase residuals PDF failed - %s\n', ex.message))
            end

            %remove figure
            close(f);
        end

        %CODE GRAPHS (with outliers)
        if plot_code == 1
            f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
            paperSize = get(f1,'PaperSize');
            set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);

            for i=1:constellations.(systems{s}).numSat

                figure(f1);
                subplot(nrows,ncols,i);
                hold on;
                title(sprintf('%c%02d',constellations.(systems{s}).sysID,constellations.(systems{s}).PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
                grid on;
                set(gca,'FontName','Verdana');
                set(gca,'FontSize',7);
                ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
                plot(RES_CODE(constellations.(systems{s}).indexes(i),:),'.b');
                hold on
                plot(find(outliers_CODE(constellations.(systems{s}).indexes(i),:) == 1),RES_CODE(constellations.(systems{s}).indexes(i), outliers_CODE(constellations.(systems{s}).indexes(i),:) == 1),'.r');
            end
            try
                print(f1 , '-dpdf', [filerootOUT '_' systems{s} '_CODE_residuals_with_OUTLIERS']);
            catch ex
                log = Logger.getInstance();
                log.addError(sprintf('saving code residuals PDF failed - %s\n', ex.message))
            end
            close(f1);
        end
    end
end

% %% GLONASS
% if constellations.GLONASS.enabled == 1
%
%
%     %PHASE GRAPHS
%     if plot_phase == 1
%         f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f,'PaperSize');
%         set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
%
%         for i=1:constellations.GLONASS.numSat
%
%             figure(f);
%             subplot(6,4,i);
%             hold on;
%             title(sprintf('R%02d',constellations.GLONASS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_PHASE(constellations.GLONASS.indexes(i),:)*1000,'.b');
%             hold on
%             plot(find(outliers_PHASE(constellations.GLONASS.indexes(i),:) == 1),RES_PHASE(constellations.GLONASS.indexes(i), outliers_PHASE(constellations.GLONASS.indexes(i),:) == 1)*1000,'.r');
%         end
%
%         print(f , '-dpdf', [filerootOUT '_GLONASS_PHASE_residuals']);
%         %remove figure
%         close(f);
%     end
%
%
%     if plot_code == 1
%         %CODE GRAPHS
%         f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f1,'PaperSize');
%         set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
%
%         for i=1:constellations.GLONASS.numSat
%             figure(f1);
%             subplot(6,4,i);
%             hold on;
%             title(sprintf('R%02d',constellations.GLONASS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_CODE(constellations.GLONASS.indexes(i),:),'.b');
%             hold on
%             plot(find(outliers_CODE(constellations.GLONASS.indexes(i),:) == 1),RES_CODE(constellations.GLONASS.indexes(i), outliers_CODE(constellations.GLONASS.indexes(i),:) == 1),'.r');
%         end
%         print(f1, '-dpdf', [filerootOUT '_GLONASS_CODE_residuals']);
%
%         close(f1);
%     end
% end
%
%
%
%
% %% Galileo
% if constellations.Galileo.enabled == 1
%
%
%     %PHASE GRAPHS
%     if plot_phase == 1
%         f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f,'PaperSize');
%         set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
%
%         for i=1:constellations.Galileo.numSat
%             figure(f);
%             subplot(6,5,i);
%             hold on;
%             title(sprintf('E%02d',constellations.Galileo.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_PHASE(constellations.Galileo.indexes(i),:)*1000,'.b');
%             hold on
%             plot(find(outliers_PHASE(constellations.Galileo.indexes(i),:) == 1),RES_PHASE(constellations.Galileo.indexes(i), outliers_PHASE(constellations.Galileo.indexes(i),:) == 1)*1000,'.r');
%         end
%
%         print(f , '-dpdf', [filerootOUT '_Galileo_PHASE_residuals']);
%         %remove figure
%         close(f);
%     end
%
%
%     %CODE GRAPHS
%     if plot_code == 1
%         f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f1,'PaperSize');
%         set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
%
%         for i=1:constellations.Galileo.numSat
%             %CODE GRAPHS
%             figure(f1);
%             subplot(6,5,i);
%             hold on;
%             title(sprintf('E%02d',constellations.Galileo.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_CODE(constellations.Galileo.indexes(i),:),'.b');
%             hold on
%             plot(find(outliers_CODE(constellations.Galileo.indexes(i),:) == 1),RES_CODE(constellations.Galileo.indexes(i), outliers_CODE(constellations.Galileo.indexes(i),:) == 1),'.r');
%         end
%         print(f1, '-dpdf', [filerootOUT '_Galileo_CODE_residuals']);
%
%         close(f1);
%     end
%
% end
%
%
%
% %% Galileo
% if constellations.BeiDou.enabled == 1
%
%
%     %PHASE GRAPHS
%     if plot_phase == 1
%         f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f,'PaperSize');
%         set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
%
%         for i=1:constellations.BeiDou.numSat
%             figure(f);
%             subplot(7,6,i);
%             hold on;
%             title(sprintf('C%02d',constellations.BeiDou.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_PHASE(constellations.BeiDou.indexes(i),:)*1000,'.b');
%             hold on
%             plot(find(outliers_PHASE(constellations.BeiDou.indexes(i),:) == 1),RES_PHASE(constellations.BeiDou.indexes(i), outliers_PHASE(constellations.BeiDou.indexes(i),:) == 1)*1000,'.r');
%         end
%
%         print(f , '-dpdf', [filerootOUT '_BeiDou_PHASE_residuals']);
%         %remove figure
%         close(f);
%     end
%
%     %CODE GRAPHS
%     if plot_code == 1
%         f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f1,'PaperSize');
%         set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
%
%         for i=1:constellations.BeiDou.numSat
%
%             figure(f1);
%             subplot(7,6,i);
%             hold on;
%             title(sprintf('C%02d',constellations.BeiDou.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_CODE(constellations.BeiDou.indexes(i),:),'.b');
%             hold on
%             plot(find(outliers_CODE(constellations.BeiDou.indexes(i),:) == 1),RES_CODE(constellations.BeiDou.indexes(i), outliers_CODE(constellations.BeiDou.indexes(i),:) == 1),'.r');
%         end
%         print(f1, '-dpdf', [filerootOUT '_BeiDou_CODE_residuals']);
%         close(f1)
%     end;
%
% end
%
%
%
% %% QZSS
% if constellations.QZSS.enabled == 1
%
%
%     %PHASE GRAPHS
%     if plot_phase == 1
%         f = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f,'PaperSize');
%         set(f,'PaperPosition',[1,1,paperSize(1)-1,paperSize(2)-1]);
%
%         for i=1:constellations.QZSS.numSat
%             %PHASE GRAPHS
%             figure(f);
%             subplot(2,2,i);
%             hold on;
%             title(sprintf('J%02d',constellations.QZSS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Phase Residual (mm)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_PHASE(constellations.QZSS.indexes(i),:)*1000,'.b');
%             hold on
%             plot(find(outliers_PHASE(constellations.QZSS.indexes(i),:) == 1),RES_PHASE(constellations.QZSS.indexes(i), outliers_PHASE(constellations.QZSS.indexes(i),:) == 1)*1000,'.r');
%         end
%
%         print(f , '-dpdf', [filerootOUT '_QZSS_PHASE_residuals']);
%         %remove figure
%         close(f);
%     end
%
%     %CODE GRAPHS
%     if plot_code == 1
%         f1 = figure('Name','goGPS processing report','NumberTitle','off','PaperOrientation','landscape','PaperUnits','centimeters','PaperType','A3','Visible','off');
%         paperSize = get(f1,'PaperSize');
%         set(f1,'PaperPosition',[2,2,paperSize(1)-2,paperSize(2)-2]);
%
%         for i=1:constellations.QZSS.numSat
%             figure(f1);
%             subplot(2,2,i);
%             hold on;
%             title(sprintf('C%02d',constellations.QZSS.PRN(i)),'FontName','Verdana','FontSize',10,'FontWeight','Bold','Color',[0 0 1]);
%             grid on;
%             set(gca,'FontName','Verdana');
%             set(gca,'FontSize',7);
%             ylabel('DD Code Residual (m)','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             xlabel('Epoch','FontName','Verdana','FontSize',6,'FontWeight','Bold');
%             plot(RES_CODE(constellations.QZSS.indexes(i),:),'.b');
%             hold on
%             plot(find(outliers_CODE(constellations.QZSS.indexes(i),:) == 1),RES_CODE(constellations.QZSS.indexes(i), outliers_CODE(constellations.QZSS.indexes(i),:) == 1),'.r');
%         end
%         print(f1, '-dpdf', [filerootOUT '_QZSS_CODE_residuals']);
%
%         close(f1);
%     end
%
% end
