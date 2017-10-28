classdef CodeTest < matlab.unittest.TestCase
    %Tester for unit testing for refactoring
    %   Detailed explanation goes here
    
    properties
%         TableHR
%         StorageHR
%         TableBR
%         StorageBR
        sg3P
        sh3P
        adjustedSol
        
    end
    
    methods (TestClassSetup)
        
        function [ testCase ] = calculateResults(testCase)
            clc
            close all
            
            
            
%             [ testCase.TableHR , testCase.StorageHR , objHR ] = mainFunction('Hydrate Ridge');
%             [ testCase.TableBR , testCase.StorageBR , objBR ] = mainFunction('Blake Ridge');
            

            
            [ testCase.sg3P , testCase.sh3P , testCase.adjustedSol ] = solubilityScript();

        end
        
    end
        
    methods (Test)
        
        function verifySg3P(testCase)
            load('Unit Testing\unitSolubility.mat');
            testCase.verifyEqual( testCase.sg3P , savedVariables.sg3P );
        end
        function verifySh3P(testCase)
            load('Unit Testing\unitSolubility.mat');
            testCase.verifyEqual( testCase.sh3P , savedVariables.sh3P );
                end
        function verifySol(testCase)
            load('Unit Testing\unitSolubility.mat');
            testCase.verifyEqual( testCase.adjustedSol , savedVariables.adjustedSol );
        end
%         function verifyStorageHR(testCase)
%             load('Unit Testing\verifyStorageHR.mat');
%             testCase.verifyEqual( unitStorageHR , testCase.StorageHR );
%         end
%         
%         function verifyTableHR(testCase)
%             load('Unit Testing\verifyTableHR.mat');
%             testCase.verifyTrue( isequaln( unitTableHR , testCase.TableHR ));
%         end
%         
%         function verifyStorageBR(testCase)
%             load('Unit Testing\verifyStorageBR.mat');
%             testCase.verifyEqual( unitStorageBR , testCase.StorageBR );
%         end
%         
%         function verifyTableBR(testCase)
%             load('Unit Testing\verifyTableBR.mat');
%             testCase.verifyTrue( isequaln( unitTableBR , testCase.TableBR ));
%         end
    end
    
end

