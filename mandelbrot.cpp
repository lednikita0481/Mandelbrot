#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <immintrin.h>

//#define NO_DRAW 1

typedef unsigned char BYTE;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHTH = 800;
const float X_MAX = 2.0;
const float Y_MAX = 2.0;
const float MAX_DISTANCE = 10.0;
const float DX = 2*X_MAX/WINDOW_WIDTH;
const float DY = 2*Y_MAX/WINDOW_HEIGHTH;

void Draw_Mandelbrot(sf::Texture* Google_Pixel, float x_mov, float y_mov, float scale);
void Draw_Mandelbrot_AVX(sf::Texture* Google_Pixel, float x_mov, float y_mov, float scale);
sf::Text *Set_Text (sf::Font &font, float x_coord, float y_coord);

int main()
{
    //printf("%d %d %f %f %f %f %f", WINDOW_WIDTH, WINDOW_HEIGHTH, X_MAX, Y_MAX, MAX_DISTANCE, DX, DY);
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHTH), "Mandelbrot!", sf::Style::Close);
    printf("%d %d %f %f %f %f %f\n", WINDOW_WIDTH, WINDOW_HEIGHTH, X_MAX, Y_MAX, MAX_DISTANCE, DX, DY);

    float x_mov = 0;
    float y_mov = 0;
    float scale = 1.0;

    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHTH);
    Draw_Mandelbrot_AVX(&texture, x_mov, y_mov, scale);
    sf::RectangleShape mandelbrot(sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHTH));
    mandelbrot.setTexture(&texture);

    sf::Clock clock; // starts the clock

    sf::Font font;
    font.loadFromFile("Disket-Mono-Bold.ttf");
    sf::Text fps_text = *Set_Text (font, 0, 0);

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event)) 
        {  
            if (event.type == sf::Event::Closed)
                window.close();
        }

        clock.restart();

        if (event.type == sf::Event::Closed)
            window.close();
        
        float x_mov_old = x_mov;
        float y_mov_old = y_mov;
        float scale_old = scale;
        
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
            x_mov -= 0.1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
            x_mov += 0.1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
            y_mov += 0.1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
            y_mov -= 0.1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::F2))
        {
            scale *= 1.25f;
            x_mov /= 1.25f;
            y_mov /= 1.25f;
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::F1))
        {
            scale /= 1.25f;
            x_mov *= 1.25f;
            y_mov *= 1.25f;
        }
        //if (x_mov_old != x_mov || y_mov_old != y_mov || scale_old != scale)
        //{   
            if (x_mov >= 2*X_MAX/scale) x_mov -= 2*X_MAX/scale;
            if (y_mov >= 2*Y_MAX/scale) y_mov -= 2*Y_MAX/scale;
            if (x_mov <= -2*X_MAX/scale) x_mov += 2*X_MAX/scale;
            if (y_mov <= -2*Y_MAX/scale) y_mov += 2*Y_MAX/scale;
            Draw_Mandelbrot_AVX(&texture, x_mov, y_mov, scale);
        //}

        sf::Time elapsed_time = clock.getElapsedTime();

        char FPS_Text[20];
        sprintf (FPS_Text, "FPS: %.2f\n", 1/elapsed_time.asSeconds());
        printf("%s\n", FPS_Text);
        fps_text.setString(FPS_Text);

        #ifndef NO_DRAW
        window.draw(mandelbrot);
        window.draw(fps_text);
        window.display();
        #endif
    }
}

void Draw_Mandelbrot(sf::Texture* Google_Pixel, float x_mov, float y_mov, float scale)
{
    #ifndef NO_DRAW
    //sf::Uint8 pixels[4*WINDOW_WIDTH*WINDOW_HEIGHTH] = {};
    sf::Image image = Google_Pixel->copyToImage();
    #endif


    float NEW_X_MAX = scale*X_MAX;
    float NEW_Y_MAX = scale*Y_MAX;
    float NEW_DX = 2*NEW_X_MAX/WINDOW_WIDTH;
    float NEW_DY = 2*NEW_Y_MAX/WINDOW_HEIGHTH;

    for (int y0_window = 0; y0_window < WINDOW_HEIGHTH; y0_window++)
    {
        float y0 = (((2*NEW_Y_MAX)) * ((float)y0_window/(float)WINDOW_HEIGHTH)) - NEW_Y_MAX + scale*y_mov;
        if (y0 > Y_MAX) y0 -= 2*Y_MAX;
        if (y0 < -Y_MAX) y0 += 2*Y_MAX;

        for (int x0_window = 0; x0_window < WINDOW_WIDTH; x0_window++)
        {
            float x0 = (((2*NEW_X_MAX))* ((float)x0_window / (float)WINDOW_WIDTH)) - NEW_X_MAX - x_mov;

            if (x0 > X_MAX) x0 -= 2*X_MAX;
            if (x0 < -X_MAX) x0 += 2*X_MAX;

            float x = x0;
            float y = y0;
            
            int cur_iter = 1;
            for (; cur_iter < 255; cur_iter++)
            {
                /*
                float x_next = (x*x) - (y*y) + x0;
                float y_next = (2*x*y) + y0;

                if ((x_next*x_next) + (y_next*y_next) > (MAX_DISTANCE*MAX_DISTANCE)) break;

                x = x_next;
                y = y_next;
                */

               float x2 = x * x, y2 = y * y, xy = x * y;
                float distance = x2 + y2;
                if (distance > MAX_DISTANCE*MAX_DISTANCE)
                    break;

                x = x2 - y2 + x0;
                y = xy + xy + y0;
            }

            //int i = 4*((WINDOW_WIDTH * y0_window) + x0_window);
            //pixels[i] = (sf::Uint8) 2*cur_iter;
            //pixels[i+1] = (sf::Uint8) 2*cur_iter;
            //pixels[i+2] = (sf::Uint8) 2*cur_iter;
            //pixels[i+3] = (sf::Uint8) 255;

            #ifndef NO_DRAW
            image.setPixel(x0_window, y0_window, sf::Color{2*cur_iter, 2*cur_iter, 2*cur_iter, 255});
            #endif
        }
    }

    #ifndef NO_DRAW
    //Google_Pixel->update(pixels);
    Google_Pixel->update(image);
    #endif
}

void Draw_Mandelbrot_AVX(sf::Texture* Google_Pixel, float x_mov, float y_mov, float scale)
{
    #ifndef NO_DRAW
    sf::Uint8 pixels[4*WINDOW_WIDTH*WINDOW_HEIGHTH] = {};
    //sf::Image image = Google_Pixel->copyToImage();
    #endif

    float NEW_X_MAX = scale*X_MAX;
    float NEW_Y_MAX = scale*Y_MAX;
    float NEW_DX = 2*NEW_X_MAX/WINDOW_WIDTH;
    float NEW_DY = 2*NEW_Y_MAX/WINDOW_HEIGHTH;

    __m256 y_max_avx = _mm256_set1_ps(Y_MAX);
    __m256 min_y_max_avx = _mm256_set1_ps(-Y_MAX);
    __m256 y_max_2 = _mm256_set1_ps(2*Y_MAX);
    __m256 x_max_avx = _mm256_set1_ps(X_MAX);
    __m256 min_x_max_avx = _mm256_set1_ps(-X_MAX);
    __m256 x_max_avx_2 = _mm256_set1_ps(2*X_MAX);
    __m256 max_distance_2 = _mm256_set1_ps(MAX_DISTANCE*MAX_DISTANCE);

    for (int y0_window = 0; y0_window < WINDOW_HEIGHTH; y0_window++)
    {
        __m256 y0_avx = _mm256_set1_ps((((2*NEW_Y_MAX)) * ((float)y0_window/(float)WINDOW_HEIGHTH)) - NEW_Y_MAX + scale*y_mov);
        

        __m256 cmp_y0_greater_avx = _mm256_cmp_ps(y0_avx, y_max_avx, _CMP_GT_OQ);
        int cmp_mask_y0_greater = _mm256_movemask_ps(cmp_y0_greater_avx);
        if (cmp_mask_y0_greater != 0)
        {
            y0_avx = _mm256_sub_ps(y0_avx, y_max_2);
        }

        __m256 cmp_y0_lower_avx = _mm256_cmp_ps(y0_avx, min_y_max_avx, _CMP_LT_OQ);
        int cmp_mask_y0_lower = _mm256_movemask_ps(cmp_y0_lower_avx);
        if (cmp_mask_y0_lower != 0)
        {
            y0_avx = _mm256_add_ps(y0_avx, y_max_2);
        }


        for (int x0_window = 0; x0_window < WINDOW_WIDTH; x0_window+=8)
        {
            __m256 x0_avx = _mm256_set1_ps((((2*NEW_X_MAX))* ((float)x0_window / (float)WINDOW_WIDTH)) - NEW_X_MAX - scale*x_mov);
            //__m256 x0_avx_offset = _mm256_setr_ps(7*DX, 6*DX, 5*DX, 4*DX, 3*DX, 2*DX, DX, 0);
            __m256 x0_avx_offset = _mm256_setr_ps(0, NEW_DX, 2*NEW_DX, 3*NEW_DX, 4*NEW_DX, 5*NEW_DX, 6*NEW_DX, 7*NEW_DX);

            x0_avx = _mm256_add_ps(x0_avx, x0_avx_offset);

            __m256 cmp_x0_greater_avx = _mm256_cmp_ps(x0_avx, x_max_avx, _CMP_GT_OQ);
            int cmp_mask_x0_greater = _mm256_movemask_ps(cmp_x0_greater_avx);
            if (cmp_mask_x0_greater != 0)
            {
                x0_avx = _mm256_sub_ps(x0_avx, x_max_avx_2);
            }

            __m256 cmp_x0_lower_avx = _mm256_cmp_ps(x0_avx, min_x_max_avx, _CMP_LT_OQ);
            int cmp_mask_x0_lower = _mm256_movemask_ps(cmp_x0_lower_avx);
            if (cmp_mask_x0_lower != 0)
            {
                x0_avx = _mm256_add_ps(x0_avx, x_max_avx_2);
            }

            __m256 x_avx = x0_avx;
            __m256 y_avx = y0_avx;

            int cur_iter = 0;
            __m256i iterations_avx = _mm256_set1_epi32(0);

            for (; cur_iter < 255; cur_iter++)
            {
                /*
                __m256 x_next_avx = _mm256_sub_ps(_mm256_mul_ps(x_avx, x_avx), _mm256_mul_ps(y_avx, y_avx));
                x_next_avx = _mm256_add_ps(x_next_avx, x0_avx);

                __m256 y_next_avx = _mm256_mul_ps(x_avx, y_avx);
                y_next_avx = _mm256_mul_ps(y_next_avx, _mm256_set1_ps(2));
                y_next_avx = _mm256_add_ps(y_next_avx, y0_avx);

                __m256 cmp_iter = _mm256_cmp_ps(_mm256_add_ps(_mm256_mul_ps(x_next_avx, x_next_avx), _mm256_mul_ps(y_next_avx, y_next_avx)), 
                                                max_distance_2, _CMP_LT_OS);
                
                int cmp_mask_iter = _mm256_movemask_ps(cmp_iter);

                if (cmp_mask_iter == 0) break;

                __m256i iter_sub = _mm256_castps_si256(cmp_iter);  // 11111111 -> -1 from mask to int (not unsigned, so we sub)

                iterations_avx = _mm256_sub_epi32 (iterations_avx, iter_sub);

                x_avx = x_next_avx;
                y_avx = y_next_avx;
                */

               __m256 x2 = _mm256_mul_ps(x_avx, x_avx);
                __m256 y2 = _mm256_mul_ps(y_avx, y_avx);
                __m256 xy = _mm256_mul_ps(x_avx, y_avx);

                __m256 dist = _mm256_add_ps(x2, y2);
                __m256 mask = _mm256_cmp_ps(dist, max_distance_2, _CMP_LT_OQ);      // FFFFFFFF (= -1) if true, 0 if false

                int res = _mm256_movemask_ps(mask);
                if (!res)  break;                                             // all distances are out of range
           
                iterations_avx = _mm256_sub_epi32 (iterations_avx, _mm256_castps_si256 (mask));   //cur_iter + 1 or cur_iter + 0

                x_avx = _mm256_add_ps(_mm256_sub_ps(x2, y2), x0_avx);
                y_avx = _mm256_add_ps(_mm256_add_ps(xy, xy), y0_avx);
            }

            int* iter_array = (int*) &iterations_avx;


            #ifndef NO_DRAW
            for (int j = 0; j < 8; j++)
            {
                int i = 4*((WINDOW_WIDTH * y0_window) + x0_window) + 4*j;
                pixels[i] = (sf::Uint8) 2*iter_array[j];
                pixels[i+1] = (sf::Uint8) 2*iter_array[j];
                pixels[i+2] = (sf::Uint8) 2*iter_array[j];
                pixels[i+3] = (sf::Uint8) 255;

                //image.setPixel(x0_window + j, y0_window, sf::Color{2*iter_array[j], 2*iter_array[j], 2*iter_array[j], 255});
            }
            #endif

        }

    }

    #ifndef NO_DRAW
    Google_Pixel->update(pixels);
    #endif
    //Google_Pixel->update(image);
}

sf::Text *Set_Text (sf::Font &font, float x_coord, float y_coord) {
    sf::Text *text = new sf::Text;          

    text->setFont(font);
    text->setFillColor(sf::Color::Yellow);
    text->setCharacterSize(30);
    text->setPosition(x_coord, y_coord);

    return text;
}